import os
from os.path import join
from functools import partial
from nipype import Function, IdentityInterface, Node, Workflow
from nipype.interfaces.io import DataSink
from nipype.interfaces.fsl.maths import MathsCommand, MeanImage

dwi_folder = "/Volumes/bronkodata_work/HCP-D/imagingcollection01/HCD0001305_V1_MR/unprocessed/Diffusion"
base_dir = "/Volumes/bronkodata_work/hcp-d-nipype"
output_dir = "/Volumes/bronkodata_work/hcp-d-nipype/output"  # Specify where you want the DataSink to store the outputs
os.makedirs(output_dir, exist_ok=True)

f_dwi = join(dwi_folder, "{subid}_{sequence}_{PEdir}.nii.gz")
subids = ["HCD0001305_V1_MR"]
sequences = ["dMRI_dir98", "dMRI_dir99"]
PEdirs = ["AP", "PA"]

dwi_imgs = [
    f_dwi.format(subid=subid, sequence=sequence, PEdir=PEdir)
    for subid in subids
    for sequence in sequences
    for PEdir in PEdirs
]


def process_images(images):
    import nibabel as nib

    from hcpy.utils import get_pedir

    pos_vols = []
    neg_vols = []

    for img_path in images:
        img = nib.load(img_path)
        if len(img.shape) < 4:
            continue
        volumes = img.shape[3]
        pedir = get_pedir(img_path)  # You'll need to define or import this function
        if "-" in pedir:
            neg_vols.append(volumes)
        else:
            pos_vols.append(volumes)

    pos_vols.sort()
    neg_vols.sort()

    return pos_vols, neg_vols


def write_volumes(pos_vols, neg_vols, path):
    from os.path import join

    pos_file_path = join(path, "pos_volumes.txt")
    neg_file_path = join(path, "neg_volumes.txt")
    pos_corr_file_path = join(path, "pos_corr_volumes.txt")
    neg_corr_file_path = join(path, "neg_corr_volumes.txt")
    with open(pos_file_path, "w") as pos_file, open(
        neg_file_path, "w"
    ) as neg_file, open(pos_corr_file_path, "w") as pos_corr_file, open(
        neg_corr_file_path, "w"
    ) as neg_corr_file:

        for vol_pos, vol_neg in zip(pos_vols, neg_vols):
            min_vol = min(vol_pos, vol_neg)
            pos_file.write(f"{min_vol} {vol_pos}\n")
            neg_file.write(f"{min_vol} {vol_neg}\n")

            if vol_pos > 0:
                pos_corr_file.write(f"{min_vol}\n")

            if vol_neg > 0:
                neg_corr_file.write(f"{min_vol}\n")

    return (
        pos_file_path,
        neg_file_path,
        pos_corr_file_path,
        neg_corr_file_path,
    )


inputnode = Node(IdentityInterface(fields=["imgs", "image_list"]), name="inputnode")
inputnode.inputs.image_list = dwi_imgs
inputnode.iterables = [("imgs", dwi_imgs)]

outputnode = Node(
    IdentityInterface(fields=["pos", "neg", "pos_corr", "neg_corr", "done"]),
    name="outputnode",
)

workflow = Workflow(name="create_series_files_workflow")
workflow.base_dir = base_dir

# Initialize DataSink
datasink = Node(DataSink(), name="datasink")
datasink.inputs.base_directory = output_dir


process_node = Node(
    Function(
        input_names=["images"],
        output_names=["pos_vols", "neg_vols"],
        function=process_images,
    ),
    name="process_images",
)
workflow.connect(inputnode, "image_list", process_node, "images")

write_volumes_node = Node(
    Function(
        input_names=[
            "pos_vols",
            "neg_vols",
            "path",
        ],
        output_names=[
            "pos_file_path",
            "neg_file_path",
            "pos_corr_file_path",
            "neg_corr_file_path",
        ],
        function=write_volumes,
    ),
    name="write_volumes",
)
write_volumes_node.inputs.path = output_dir

workflow.connect(
    [
        (
            process_node,
            write_volumes_node,
            [("pos_vols", "pos_vols"), ("neg_vols", "neg_vols")],
        ),
        (
            write_volumes_node,
            datasink,
            [
                ("pos_file_path", "@pos_volumes"),
                ("neg_file_path", "@neg_volumes"),
                ("pos_corr_file_path", "@pos_corr_volumes"),
                ("neg_corr_file_path", "@neg_corr_volumes"),
            ],
        ),
    ]
)

# Node to calculate mean along X, Y, and Z (spatial mean), then T (temporal mean)
calc_mean = Node(
    MeanImage(args="-Xmean -Ymean -Zmean -Tmean"),
    name="calc_mean",
)
workflow.connect(inputnode, "imgs", calc_mean, "in_file")
workflow.connect(calc_mean, "out_file", datasink, "rawdata.@mean_volumes")


def extract_b0_volumes(nii_path, b0maxbval, b0dist):
    import os
    import nibabel as nib
    from nipype.interfaces.fsl import ExtractROI, Merge
    import subprocess
    from nilearn.image import index_img
    from hcpy.utils import extract_b0_indices, get_pedir, get_readout_time

    # Define or import your extract_b0_indices, get_pedir, and get_readout_time functions here
    bval_path = nii_path.replace(".nii.gz", ".bval")

    # Get b0 indices
    b0_indices, log = extract_b0_indices(bval_path, b0maxbval, b0dist=b0dist)
    b0_dwi = index_img(nii_path, b0_indices).to_filename(
        nii_path.replace(".nii.gz", "_b0.nii.gz")
    )

    return nii_path.replace(".nii.gz", "_b0.nii.gz")


extract_b0_node = Node(
    Function(
        input_names=[
            "nii_path",
            "b0maxbval",
            "b0dist",
            "logdir",
        ],
        output_names=["b0_image"],
        function=extract_b0_volumes,
    ),
    name="extract_b0_volumes",
)
extract_b0_node.inputs.b0maxbval = 50
extract_b0_node.inputs.b0dist = 45

workflow.connect(inputnode, "imgs", extract_b0_node, "nii_path")
workflow.connect(extract_b0_node, "b0_image", datasink, "rawdata.@b0_volumes")

workflow.write_graph(graph2use="flat", dotfilename="workflow_graph.dot")
workflow.run()
