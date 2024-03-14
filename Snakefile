# basic_preproc_norm_intensity.snakefile
import os
import json
from os.path import join, basename, dirname
import hcpy
from hcpy.utils import *


# Define global variables and configurations
configfile: "config.yaml"


topup_conf = hcpy.config.b02b0

f_dwi = join(
    config["rawdir"],
    "{subid}",
    "unprocessed",
    "Diffusion",
    "{subid}_{sequence}_{PEdir}.nii.gz",
)
f_bval = join(
    config["rawdir"],
    "{subid}",
    "unprocessed",
    "Diffusion",
    "{subid}_{sequence}_{PEdir}.bval",
)
f_bvec = join(
    config["rawdir"],
    "{subid}",
    "unprocessed",
    "Diffusion",
    "{subid}_{sequence}_{PEdir}.bvec",
)

subdir = join(config["outdir"], "{subid}")
qcdir = join(subdir, "Diffusion", "qc")
topupdir = join(subdir, "Diffusion", "topup")
eddydir = join(subdir, "Diffusion", "eddy")
datadir = join(subdir, "Diffusion", "data")
rawdatadir = join(subdir, "Diffusion", "rawdata")
logdir = join(subdir, "Diffusion", "log")
outdir_T1w = join(subdir, "T1w", "Diffusion")

logfiles = [
    join(logdir, "combine_data.done"),
    join(logdir, "get_rotated_bvecs.done"),
    join(logdir, "average_bvecs.done"),
    join(logdir, "create_fov_mask.done"),
    join(logdir, "mask_out_data.done"),
    join(logdir, "epi_reg_dof.done"),
    join(logdir, "remove_bias_field.done"),
    join(logdir, "bbregister.done"),
    join(logdir, "tkregister2.done"),
    join(logdir, "apply_warp3_remove_bias.done"),
    join(logdir, "register_diff_to_T1w.done"),
]

f_b0 = join(rawdatadir, "{subid}_{sequence}_{PEdir}_b0.nii.gz")

f_dwi_res = join(rawdatadir, "{subid}_{sequence}_{PEdir}_res.nii.gz")

f_checkeven = join(config["logdir"], "{subid}", "{subid}_{sequence}_{PEdir}.checkeven")
checkfiles = [f_checkeven]

# QC

f_qc_topup = join(qcdir, "{subid}_topup.svg")
f_qc_bvecs = join(qcdir, "{subid}_bvecs.html")


# topup
f_b0_log = join(topupdir, "extractedb0.txt")
# f_b0_index = join(topupdir, "b0_index.txt")
f_b0_posneg = join(topupdir, "Pos_Neg_b0.nii.gz")
f_b0_pos = join(topupdir, "Pos_b0.nii.gz")
f_b0_neg = join(topupdir, "Neg_b0.nii.gz")
f_acqparams_topup = join(topupdir, "acqparams.txt")
f_topup_field = join(topupdir, "topup_Pos_Neg_b0_field.nii.gz")
f_b0_hifi = join(topupdir, "hifib0.nii.gz")

# eddy
f_index = join(eddydir, "index.txt")
f_series_index = join(eddydir, "series_index.txt")
f_img_posneg = join(eddydir, "Pos_Neg.nii.gz")
f_bvals_posneg = join(eddydir, "Pos_Neg.bvals")
f_bvecs_posneg = join(eddydir, "Pos_Neg.bvecs")
f_acqparams_eddy = join(eddydir, "acqparams.txt")
f_eddy = join(eddydir, "eddy_unwarped_images.nii.gz")


sequences = config["ImagingSequences"]
PEdirs = config["PhaseEncodingDirections"]
subids = [
    "HCD0001305_V1_MR",
]


outfiles = [
    join(outdir_T1w, "nodif_brain_mask_new.nii.gz"),
    join(outdir_T1w, "data_masked_thresh.nii.gz"),
    join(outdir_T1w, "bvals"),
    join(outdir_T1w, "bvecs"),
]
qcfiles = [
    join(logdir, "qc1.done"),
    # join(logdir, "qc2.done"),
    join(logdir, "qc3.done"),
    join(logdir, "qc4.done"),
]
logfiles = [
    join(logdir, f)
    for f in [
        "tkregister2.done",
        "bbregister.done",
        "get_transformation_matrices.done",
        "apply_warp3.done",
        "remove_bias.done",
        "diff_res_structural_space.done",
        "diff_res_mask.done",
        "dilate_mask.done",
        "rotate_bvecs_to_struct.done",
        "register_diff_to_T1w.done",
        "adjust_fov_mask.done",
        "fov_mask_data.done",
        "threshold_data.done",
        "new_nodif_brain_mask.done",
        "calculate_coverage.done",
        "create_mif.done",
        "tissue_segmentation.done",
        "bias_correction.done",
    ]
]


rule all:
    """
    Generates all required outputs.
    """
    input:
        expand(
            qcfiles + logfiles,
            subid=subids,
            sequence=sequences,
            PEdir=PEdirs,
        ),


rule create_series_files:
    input:
        imgs=expand([f_dwi], subid="{subid}", sequence=sequences, PEdir=PEdirs),
    output:
        pos=join(eddydir, "Pos_SeriesVolNum.txt"),
        neg=join(eddydir, "Neg_SeriesVolNum.txt"),
        pos_corr=join(subdir, "Diffusion", "rawdata", "Pos_SeriesCorrespVolNum.txt"),
        neg_corr=join(subdir, "Diffusion", "rawdata", "Neg_SeriesCorrespVolNum.txt"),
    run:
        import json
        import nibabel as nib

        # Initialize lists to store volume numbers
        pos_vols = []
        neg_vols = []

        # Process each image to categorize them and count volumes
        for img in input.imgs:
            volumes = nib.load(img).shape[3]
            pedir = get_pedir(img)
            if "-" in pedir:
                neg_vols.append(volumes)
            else:
                pos_vols.append(volumes)

                # Ensure volumes are sorted to match pairs correctly
        pos_vols.sort()
        neg_vols.sort()

        with open(output.pos, "w") as pos_file, open(output.neg, "w") as neg_file, open(
            output.pos_corr, "w"
        ) as pos_corr_file, open(output.neg_corr, "w") as neg_corr_file:
            for vol_pos, vol_neg in zip(pos_vols, neg_vols):
                min_vol = min(vol_pos, vol_neg)
                pos_file.write(f"{min_vol} {vol_pos}\n")
                neg_file.write(f"{min_vol} {vol_neg}\n")

                if vol_pos > 0:
                    pos_corr_file.write(f"{min_vol}\n")

                if vol_neg > 0:
                    neg_corr_file.write(f"{min_vol}\n")


rule verify_input_pairs:
    """
    Verifies that input datasets are provided in pairs.
    """
    input:
        pos=expand(
            f_dwi,
            rawdir=config["rawdir"],
            subid="{subid}",
            sequence=sequences,
            PEdir=PEdirs[0],
        ),
        neg=expand(
            f_dwi,
            rawdir=config["rawdir"],
            subid="{subid}",
            sequence=sequences,
            PEdir=PEdirs[1],
        ),
    run:
        pos_count = len(input.pos)
        neg_count = len(input.neg)
        if pos_count != neg_count:
            raise ValueError(
                "Wrong number of input datasets! Make sure that you provide pairs of input filenames."
            )
        else:
            print("Input datasets are provided in pairs.")


rule check_even_slices:
    """
    Load f_dwi with nibabel and check if the number of slices is even.
    """
    input:
        dwi=f_dwi,
    output:
        checkeven=f_checkeven,
    run:
        import nibabel as nib

        is_even_slices = nib.load(input.dwi).shape[2] % 2 == 0
        if basename(topup_conf) == "b02b0.cnf" and not is_even_slices:
            raise ValueError(
                "Input images have an odd number of slices. This is incompatible with the default topup configuration file. Either supply a topup configuration file that doesn't use subsampling (e.g., FSL's 'b02b0_1.cnf') using the --topup-config-file=<file> flag (recommended), or instruct the HCP pipelines to remove a slice using the --ensure-even-slices flag (legacy option). Note that the legacy option is incompatible with slice-to-volume correction in FSL's eddy."
            )
        with open(output.checkeven, "w") as f:
            f.write(str(is_even_slices))


################################
# basic_preproc_norm_intensity #
################################
rule calculate_mean:
    """
    Calculates the mean of input DWI images.
    """
    input:
        raw=f_dwi,
    output:
        mean=join(rawdatadir, "{subid}_{sequence}_{PEdir}_mean.nii.gz"),
    shell:
        """
        fslmaths {input.raw} -Xmean -Ymean -Zmean {output.mean}
        fslmaths {output.mean} -Tmean {output.mean}
        """


rule meants:
    """
    Computes the mean timeseries.
    """
    input:
        mean=rules.calculate_mean.output.mean,
    output:
        meants=join(rawdatadir, "{subid}_{sequence}_{PEdir}_meants.txt"),
    shell:
        """
        echo $(fslmeants -i {input.mean}) >> {output.meants}
        """


rule extract_b0_volumes:
    """
    Extracts b0 volumes from the DWI data.
    """
    input:
        nii=f_dwi,
        bval=f_bval,
    output:
        b0=f_b0,
    params:
        b0maxbval=config["b0maxbval"],
        b0dist=config["b0dist"],
    run:
        import os
        import nibabel as nib

        # Get b0 indices
        indices, log = extract_b0_indices(
            input.bval, params.b0maxbval, b0dist=params.b0dist
        )

        indcount = 1  # Initialize index counter

        # Prepare FSL command to extract b0 volumes
        for idx in indices:
            output_file = os.path.splitext(output.b0)[0] + f"_{idx}.nii.gz"
            shell(f"fslroi {input.nii} {output_file} {idx} 1")
            indcount += 1

        PEdir = get_pedir(input.nii)
        ro_time = get_readout_time(input.nii, PEdir)

        # Merge the extracted b0 volumes
        shell(f"fslmerge -t {output.b0} {os.path.splitext(output.b0)[0]}_*.nii.gz")
        shell(f"rm {os.path.splitext(output.b0)[0]}_*.nii.gz")


rule rescale_series_conditionally:
    """
    The rule rescales each image in the DWI series based on the scaling factor of the first image.
    """
    input:
        dwi=f_dwi,
        json=f_dwi.replace(".nii.gz", ".json"),
        scaleS=rules.meants.output.meants,
    output:
        rescaled_dwi=f_dwi_res,
        rescaled_json=f_dwi_res.replace(".nii.gz", ".json"),
    run:
        import subprocess

        # Check if the current series and PEdir are the first ones in their respective lists
        rescale = (
            input.scaleS.replace("{rawdir}", config["rawdir"])
            .replace("{subid}", wildcards.subid)
            .replace("{PEdir}", PEdirs[0])
            .replace("{sequence}", sequences[0])
        )
        cmd = f"fslmaths {input.dwi} -mul $(cat {rescale}) -div $(cat {input.scaleS}) {output.rescaled_dwi}"
        print(cmd)
        shell(cmd)
        shell(f"cp {input.json} {output.rescaled_json}")


##########################
# basic_preproc_sequence #
##########################


rule merge_data:
    """
    Merges rescaled DWI images and related files based on phase encoding directions.

    This rule processes input DWI images, b0 images, b-values, b-vectors, acquisition parameters,
    and indices for both positive and negative phase encoding directions. It outputs merged images
    for both positive and negative directions, as well as combined images and parameters for further
    processing steps such as topup and eddy correction.

    Inputs:
        imgs: Rescaled DWI images.
        b0_imgs: Extracted b0 images.
        bvals: b-value files.
        bvecs: b-vector files.
        acqparams: Acquisition parameter files.
        b0_indices: b0 indices files.
        series_indices: Series index files.

    Outputs:
        img_pos: Merged positive phase encoding direction images.
        img_neg: Merged negative phase encoding direction images.
        b0_posneg: Combined b0 images for topup.
        b0_pos: Merged b0 images with positive phase encoding.
        b0_neg: Merged b0 images with negative phase encoding.
        img_posneg: Combined DWI images for eddy.
        bval_pos: b-values for positive phase encoding direction.
        bval_neg: b-values for negative phase encoding direction.
        bvals_posneg: Combined b-values.
        bvecs_posneg: Combined b-vectors.
    """
    input:
        imgs=expand(
            [f_dwi_res],
            subid="{subid}",
            sequence=sorted(sequences),
            PEdir=sorted(PEdirs),
        ),
        b0_imgs=expand(
            [f_b0],
            subid="{subid}",
            sequence=sorted(sequences),
            PEdir=sorted(PEdirs),
        ),
        bvals=expand(
            [f_bval], subid="{subid}", sequence=sorted(sequences), PEdir=sorted(PEdirs)
        ),
        bvecs=expand(
            [f_bval.replace(".bval", ".bvec")],
            subid="{subid}",
            sequence=sorted(sequences),
            PEdir=sorted(PEdirs),
        ),
    params:
        b0maxbval=config["b0maxbval"],
        b0dist=config["b0dist"],
    output:
        # topup
        topup_acqparams=f_acqparams_topup,
        b0_log=f_b0_log,
        b0_posneg=f_b0_posneg,
        b0_pos=f_b0_pos,
        b0_neg=f_b0_neg,
        # eddy
        eddy_acqparams=f_acqparams_eddy,
        img_posneg=f_img_posneg,
        bval_pos=join(eddydir, "Pos.bval"),
        bval_neg=join(eddydir, "Neg.bval"),
        bvec_pos=join(eddydir, "Pos.bvec"),
        bvec_neg=join(eddydir, "Neg.bvec"),
        bvals_posneg=f_bvals_posneg,
        bvecs_posneg=f_bvecs_posneg,
        series_index=f_series_index,
        index=f_index,
    run:
        import numpy as np
        import nibabel as nib
        from nilearn.image import concat_imgs

        np.set_printoptions(precision=6)

        acqparams = ""
        pos_imgs = []
        neg_imgs = []
        pos_b0s = []
        neg_b0s = []
        series_index = np.array([])
        pos_bvals = np.array([])
        neg_bvals = np.array([])
        pos_bvecs = np.empty((3, 0))
        neg_bvecs = np.empty((3, 0))
        neg_acqparams = []
        pos_acqparams = []
        neg_log = []
        pos_log = []

        ind = 1
        b0_ind = 1

        print("Number of images:", len(input.imgs))
        print("Number of b0 images:", len(input.b0_imgs))
        print("Number of b-values:", len(input.bvals))
        print("Number of b-vectors:", len(input.bvecs))

        for (
            img,
            b0_img,
            bval,
            bvec,
        ) in zip(
            input.imgs,
            input.b0_imgs,
            input.bvals,
            input.bvecs,
        ):
            series_index = np.append(
                series_index, np.repeat(ind, nib.load(img).shape[3])
            )

            acqparams, b0_volumes, log = get_b0_acqparms(
                img, bval, b0maxbval=params.b0maxbval, b0dist=params.b0dist
            )

            bval = np.loadtxt(bval)

            pedir = get_pedir(img)
            ro_time = get_readout_time(img, pedir)
            # acqparam = get_acqparams(img)
            print(f"pedir: {pedir}")
            if "-" in pedir:
                neg_imgs.append(img)
                neg_b0s.append(b0_img)
                neg_bvals = np.append(neg_bvals, bval, axis=0)
                neg_bvecs = np.append(neg_bvecs, np.loadtxt(bvec), axis=1)
                neg_log.append(log)
                neg_acqparams.extend(acqparams)
            else:
                pos_imgs.append(img)
                pos_b0s.append(b0_img)
                pos_bvals = np.append(pos_bvals, bval, axis=0)
                pos_bvecs = np.append(pos_bvecs, np.loadtxt(bvec), axis=1)
                pos_log.append(log)
                pos_acqparams.extend(acqparams)

            ind += 1

        with open(output.b0_log, "w") as f:
            f.write("\n".join(pos_log + neg_log))

        with open(output.topup_acqparams, "w") as f:
            f.write("\n".join(pos_acqparams + neg_acqparams))
        with open(output.eddy_acqparams, "w") as f:
            f.write("\n".join(pos_acqparams + neg_acqparams))

        np.savetxt(output.series_index, series_index, fmt="%i")
        print()
        print(f"neg_imgs: {neg_b0s}")
        print(f"pos_imgs: {pos_b0s}")
        print()

        shell(f"fslmerge -t {output.b0_pos} {' '.join(pos_b0s)}")
        shell(f"fslmerge -t {output.b0_neg} {' '.join(neg_b0s)}")
        shell(f"fslmerge -t {output.b0_posneg} {output.b0_pos} {output.b0_neg}")

        shell(f"fslmerge -t {output.img_posneg} {' '.join(pos_imgs + neg_imgs)}")

        np.savetxt(
            output.bval_pos, pos_bvals.reshape(-1)[None], fmt="%i", delimiter=" "
        )
        np.savetxt(
            output.bval_neg, neg_bvals.reshape(-1)[None], fmt="%i", delimiter=" "
        )
        np.savetxt(output.bvec_pos, pos_bvecs, fmt="%.6f", delimiter=" ")
        np.savetxt(output.bvec_neg, neg_bvecs, fmt="%.6f", delimiter=" ")

        bvals = np.append(pos_bvals, neg_bvals, axis=0)
        bvecs = np.append(pos_bvecs, neg_bvecs, axis=1)

        # Create index file
        index = get_eddy_index(bvals, b0maxbval=params.b0maxbval, b0dist=params.b0dist)
        np.savetxt(output.index, index, fmt="%i")

        # Save bvals and bvecs
        np.savetxt(
            output.bvals_posneg, bvals.reshape(-1)[None], fmt="%i", delimiter=" "
        )
        np.savetxt(output.bvecs_posneg, bvecs, fmt="%.6f", delimiter=" ")


rule qc_bvecs_merged:
    input:
        bvecs=f_bvecs_posneg,
        bvals=f_bvals_posneg,
    output:
        fig=f_qc_bvecs,
    run:
        from hcpy.plotting import plot_bvecs

        plot_bvecs(input.bvecs, bvals=input.bvals, fout=output.fig)


rule qc1:
    """
    Generates QC plots for DWI data.
    """
    input:
        dwis=expand(f_dwi, subid="{subid}", sequence=sequences, PEdir=PEdirs),
        rescaled_dwis=expand(
            f_dwi_res, subid="{subid}", sequence=sequences, PEdir=PEdirs
        ),
        b0s=expand(f_b0, subid="{subid}", sequence=sequences, PEdir=PEdirs),
    output:
        fig=join(qcdir, "basic_preproc_norm_intensity.svg"),
        log=join(logdir, "qc1.done"),
    run:
        from nilearn import plotting
        import matplotlib.pyplot as plt
        import matplotlib
        from nilearn.image import concat_imgs, index_img

        matplotlib.use("agg")
        fig, axs = plt.subplots(nrows=len(input.dwis), ncols=3, figsize=(16, 9))

        i = 0
        for dwi, res_dwi, b0 in zip(input.dwis, input.rescaled_dwis, input.b0s):
            # Use a non-interactive backend
            plotting.plot_stat_map(index_img(dwi, 0), axes=axs[i, 0])
            plotting.plot_stat_map(index_img(res_dwi, 0), axes=axs[i, 1])
            plotting.plot_stat_map(index_img(b0, 0), axes=axs[i, 2])
            i += 1

        fig.savefig(output.fig)
        shell("touch {output.log}")


################################
# run topup                    #
################################
rule topup:
    """
    Runs topup to estimate susceptibility-induced distortions.
    """
    input:
        b0_posneg=f_b0_posneg,
        acqparams=f_acqparams_topup,
        topup_conf=topup_conf,
    output:
        field=f_topup_field,
    shell:
        """
        topupdir=$(dirname {output.field})
        topup --imain={input.b0_posneg} --datain={input.acqparams} --config={input.topup_conf} --out=$topupdir/topup_Pos_Neg_b0 -v --fout={output.field}
        """


rule get_hifi_b0:
    input:
        field=rules.topup.output.field,
        b0_pos=f_b0_pos,
        b0_neg=f_b0_neg,
        acqparams=f_acqparams_topup,
    output:
        b0_hifi=f_b0_hifi,
    shell:
        """
        topupdir=$(dirname {input.field})
        fslroi {input.b0_pos} {input.b0_pos}01 0 1
        fslroi {input.b0_neg} {input.b0_neg}01 0 1
        applytopup --imain={input.b0_pos}01,{input.b0_neg}01 --topup=$topupdir/topup_Pos_Neg_b0 --datain={input.acqparams} --inindex=1,$(($(fslval {input.b0_pos} dim4)+ 1)) --out={output.b0_hifi}
        """


rule bet_hifi_b0:
    input:
        b0=rules.get_hifi_b0.output.b0_hifi,
    output:
        nodif_brain=join(topupdir, "nodif_brain.nii.gz"),
    shell:
        """
        bet {input.b0} {output.nodif_brain} -m -f 0.3
        """


rule qc2:
    input:
        nodif_brain=rule.bet_hifi_b0.output.nodif_brain,
        nodif_brain_mask=rule.bet_hifi_b0.output.nodif_brain.replace(
            ".nii.gz", "_mask.nii.gz"
        ),
    output:
        fig=f_qc_topup,
    run:
        from nilearn import plotting
        import matplotlib.pyplot as plt
        import matplotlib
        from nilearn.image import concat_imgs, index_img

        # Use a non-interactive backend
        matplotlib.use("agg")
        fig, axs = plt.subplots(nrows=2, figsize=(16, 9))
        plotting.plot_anat(input.nodif_brain, axes=axs[0])
        plotting.plot_roi(input.nodif_brain_mask, input.nodif_brain, axes=axs[1])
        fig.savefig(output.fig)


################################
# run eddy                     #
################################


rule eddy:
    """
    Correct for eddy current-induced distortions and subject movements
    """
    input:
        imain=rules.merge_data.output.img_posneg,
        field=rules.topup.output.field,
        mask=rules.bet_hifi_b0.output.nodif_brain,
        index=f_series_index,
        acqparams=f_acqparams_eddy,
        bvecs=f_bvecs_posneg,
        bvals=f_bvals_posneg,
    params:
        fwhm=0,
    output:
        eddy=f_eddy,
    run:
        topup = input.field.replace("_field.nii.gz", "")
        shell(
            "eddy"
            + " diffusion"
            + " --cnr_maps"
            + f" --imain={input.imain}"
            + f" --mask={input.mask}"
            + f" --index={input.index}"
            + f" --acqp={input.acqparams}"
            + f" --bvecs={input.bvecs}"
            + f" --bvals={input.bvals}"
            + f" --fwhm={params.fwhm}"
            + f" --topup={topup}"
            + f" --out={output.eddy.replace('.nii.gz', '')}"
            + f" --nthr={config['nthreads']}"
            + " --verbose"
        )


################################
# Post-Eddy                    #
################################

dof = 6
combine_data = False
select_best_b0 = False

# TODO: Add option for Gradient Nonlinearity Distortion Correction


# eddy_postproc.sh
rule qc3:
    """Eddy qc report"""
    input:
        eddy=f_eddy,
        index=f_index,
        acqparams=f_acqparams_eddy,
        mask=join(topupdir, "nodif_brain_mask.nii.gz"),
        bvals=f_bvals_posneg,
        bvecs=join(eddydir, "eddy_unwarped_images.eddy_rotated_bvecs"),
        field=f_topup_field,
    output:
        report=directory(join(qcdir, "eddy")),
        log=join(logdir, "qc3.done"),
    run:
        shell("rm -rf {output.report}")
        cmd = (
            f"eddy_quad {input.eddy.replace('.nii.gz', '')}"
            + f" -idx {input.index}"
            + f" -par {input.acqparams}"
            + f" -m {input.mask}"
            + f" -b {input.bvals}"
            + f" -g {input.bvecs}"
            + f" -f {input.field}"
            + f" -o {output.report}"
            + " -v"
        )
        shell(cmd)
        shell("touch {output.log}")


rule combine_data:
    """
    Combine data across diffusion directions with opposing phase-encoding directions. This rule is only run if the `combine_data` flag is set to True.
    """
    input:
        pos_bvals=rules.merge_data.output.bval_pos,
        neg_bvals=rules.merge_data.output.bval_neg,
        pos_bvecs=rules.merge_data.output.bvec_pos,
        neg_bvecs=rules.merge_data.output.bvec_neg,
        unwarped_imgs=join(eddydir, "eddy_unwarped_images.nii.gz"),
        neg_series_vol_num=rules.create_series_files.output.neg,
        pos_series_vol_num=rules.create_series_files.output.pos,
    output:
        data=protected(join(datadir, "data.nii.gz")),
        log=join(logdir, "combine_data.done"),
    params:
        data=directory(datadir),
        SelectBestB0=0,
        CombineDataFlag=1,
        unwarped_pos=join(eddydir, "eddy_unwarped_Pos.nii.gz"),
        unwarped_neg=join(eddydir, "eddy_unwarped_Neg.nii.gz"),
    run:
        import os

        EddyJacFlag = "JacobianResampling"
        PosVols = len(np.loadtxt(input.pos_bvals))
        NegVols = len(np.loadtxt(input.neg_bvals))

        os.makedirs(params.data, exist_ok=True)
        shell(
            f"fslroi {input.unwarped_imgs} {params.unwarped_pos} {params.SelectBestB0} {PosVols}"
        )
        shell(
            f"fslroi {input.unwarped_imgs} {params.unwarped_neg} {PosVols + params.SelectBestB0} {NegVols}"
        )
        combine_cmd = (
            f"eddy_combine {params.unwarped_pos} {input.pos_bvals} {input.pos_bvecs} {input.pos_series_vol_num}"
            + f" {params.unwarped_neg} {input.neg_bvals} {input.neg_bvecs} {input.neg_series_vol_num}"
            + f" {params.data}"
            + f" {params.CombineDataFlag}"
        )
        shell(combine_cmd)
        shell("imrm {params.unwarped_pos} {params.unwarped_neg}")

        shell(f"mv {params.data}/bvecs {params.data}/bvecs_noRot")
        shell(f"mv {params.data}/bvals {params.data}/bvals_noRot")

        shell(f"echo {PosVols}, {NegVols} >> {output.log}")


rule get_rotated_bvecs:
    input:
        pos_bvals=rules.merge_data.output.bval_pos,
        neg_bvals=rules.merge_data.output.bval_neg,
        bvecs=join(eddydir, "eddy_unwarped_images.eddy_rotated_bvecs"),
    output:
        pos_bvecs=join(eddydir, "Pos_rotated.bvec"),
        neg_bvecs=join(eddydir, "Neg_rotated.bvec"),
        log=join(logdir, "get_rotated_bvecs.done"),
    params:
        SelectBestB0=0,
        eddydir=directory(eddydir),
    run:
        import numpy as np

        PosVols = len(np.loadtxt(input.pos_bvals))
        NegVols = len(np.loadtxt(input.neg_bvals))
        bvecs = np.loadtxt(input.bvecs)
        pos_bvecs = bvecs[:, : PosVols + params.SelectBestB0]
        neg_bvecs = bvecs[:, PosVols + params.SelectBestB0 :]
        np.savetxt(output.pos_bvecs, pos_bvecs, fmt="%.6f", delimiter=" ")
        np.savetxt(output.neg_bvecs, neg_bvecs, fmt="%.6f", delimiter=" ")

        shell("touch {output.log}")


rule average_rotated_bvecs:
    input:
        bval_pos=rules.merge_data.output.bval_pos,
        bval_neg=rules.merge_data.output.bval_neg,
        bvec_pos_rot=rules.get_rotated_bvecs.output.pos_bvecs,
        bvec_neg_rot=rules.get_rotated_bvecs.output.neg_bvecs,
        pos_series_vol_num=rules.create_series_files.output.pos,
        neg_series_vol_num=rules.create_series_files.output.neg,
    params:
        eddydir=directory(eddydir),
        datadir=directory(datadir),
        CombineDataFlag=1,
    output:
        log=join(logdir, "average_bvecs.done"),
        av_bval=join(datadir, "bvals"),
        av_bvec=join(datadir, "bvecs"),
        av_indices=join(datadir, "avg_data.idxs"),
    run:
        from hcpy import average_bvecs

        average_bvecs.main(
            input.bval_pos,
            input.bvec_pos_rot,
            input.bval_neg,
            input.bvec_neg_rot,
            join(params.datadir, "avg_data.bval"),
            join(params.datadir, "avg_data.bvec"),
            indicesoutfile=output.av_indices,
            only_combine=params.CombineDataFlag == 1,
            overlap1file=input.pos_series_vol_num,
            overlap2file=input.neg_series_vol_num,
        )
        shell(f"mv {params.datadir}/avg_data.bval {output.av_bval}")
        shell(f"mv {params.datadir}/avg_data.bvec {output.av_bvec}")
        shell(f"rm -f {params.datadir}/avg_data.bval {params.datadir}/avg_data.bvec")
        shell(f"echo $(date) > {output.log}")


rule create_fov_mask:
    input:
        data=rules.combine_data.output.data,
    params:
        eddydir=directory(eddydir),
        datadir=directory(datadir),
    output:
        cnr_maps=join(datadir, "cnr_maps.nii.gz"),
        log=join(logdir, "create_fov_mask.done"),
        fov_mask=join(datadir, "fov_mask.nii.gz"),
    run:
        from os.path import join

        shell(
            f"imcp {join(params.eddydir, 'eddy_unwarped_images.eddy_cnr_maps.nii.gz')} {output.cnr_maps}"
        )

        shell(f"fslmaths {input.data} -abs -Tmin -bin -fillh {output.fov_mask}")
        shell(f"echo $(date) > {output.log}")


# TODO: Add rule for Gradient Nonlinearity Distortion Correction (if [ ! $GdCoeffs = "NONE" ])


rule mask_out_data:
    input:
        fov_mask=rules.create_fov_mask.output.fov_mask,
        data=rules.combine_data.output.data,
        cnr_maps=rules.create_fov_mask.output.cnr_maps,
    output:
        data_masked=join(datadir, "data_masked.nii.gz"),
        cnr_maps_masked=join(datadir, "cnr_maps_masked.nii.gz"),
        nodif=join(datadir, "nodif.nii.gz"),
        nodif_brain=join(datadir, "nodif_brain.nii.gz"),
        log=join(logdir, "mask_out_data.done"),
    params:
        logdir=directory(logdir),
        datadir=directory(datadir),
    run:
        import datetime

        # from hcpy.workflows.Diffusion.PostEddy import mask_and_clean
        # start1 = datetime.datetime.now()
        # wf = mask_and_clean.setup_workflow(input, params, output)
        # wf.run("MultiProc", plugin_args={"n_procs": config["nthreads"]})
        # stop1 = datetime.datetime.now()
        # runtime = stop1 - start1
        # print(f"Runtime: {runtime}")
        # TODO: RENAME FILES AFTER WORKFLOW

        start = datetime.datetime.now()
        # mask out any data outside the field of view
        shell(f"fslmaths {input.data} -mas {input.fov_mask} {output.data_masked}")
        shell(
            f"fslmaths {input.cnr_maps} -mas {input.fov_mask} {output.cnr_maps_masked}"
        )
        # Remove negative intensity values (from eddy) from final data
        shell(f"fslmaths {output.data_masked} -thr 0 {output.data_masked}")
        shell(f"fslroi {output.data_masked} {output.nodif} 0 1")
        shell(f"bet {output.nodif} {output.nodif_brain} -m -f 0.1")
        stop = datetime.datetime.now()
        print(f"Runtime: {stop - start}")


        shell(f"echo $(date) > {output.log}")


################################
# DiffusionToStructural        #
################################
dir_T1w = join(config["t1dir"], "{subid}", "T1w")
f_T1w = join(dir_T1w, "T1w_acpc_dc.nii.gz")
f_T1w_restore = join(dir_T1w, "T1w_acpc_dc_restore.nii.gz")
f_T1w_restore_brain = join(dir_T1w, "T1w_acpc_dc_restore_brain.nii.gz")

f_BiasField = join(dir_T1w, "BiasField_acpc_dc.nii.gz")
f_FreeSurferBrainMask = join(dir_T1w, "brainmask_fs.nii.gz")

f_RegOutput = join(subdir, "Diffusion", "reg", "RegOutput")
f_QAImage = join(subdir, "Diffusion", "reg", "T1wMulEPI")
# DiffRes=$(${FSLDIR}/bin/fslval ${outdir}/data/data pixdim1)
# DiffRes=$(printf "%0.2f" ${DiffRes})


rule epi_reg_dof:
    """b0 FLIRT BBR and bbregister to T1w"""
    input:
        t1=f_T1w,
        t1_brain=f_T1w_restore_brain,
        epi=rules.mask_out_data.output.nodif,
    output:
        mat=join(outdir_T1w, "nodif2T1w_initII.mat"),
        init_mat=join(outdir_T1w, "nodif2T1w_initII_init.mat"),
        log=join(logdir, "epi_reg_dof.done"),
    params:
        dof=dof,
        basepath=join(outdir_T1w, "nodif2T1w_initII"),
    group:
        "DiffusionToStructural"
    run:
        # TODO: write out EPI_REG_DOF into workflow
        shell(
            f"$HCPPIPEDIR/global/scripts/epi_reg_dof --dof={params.dof} --epi={input.epi} --t1={input.t1} --t1brain={input.t1_brain} --out={params.basepath}"
        )

        shell(f"echo $(date) > {output.log}")


rule apply_warp1:
    input:
        epi=rules.mask_out_data.output.nodif,
        premat=rules.epi_reg_dof.output.init_mat,
        t1=f_T1w,
    output:
        epi2t1=join(outdir_T1w, "nodif2T1w_init.nii.gz"),
        log=join(logdir, "apply_warp.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(
            f"applywarp -i {input.epi} -o {output.epi2t1} -r {input.t1} --rel --interp=spline --premat={input.premat}"
        )


rule apply_warp2:
    input:
        epi=rules.mask_out_data.output.nodif,
        premat=rules.epi_reg_dof.output.mat,
        t1=f_T1w,
    output:
        epi2t1=join(outdir_T1w, "nodif2T1w_initII.nii.gz"),
        log=join(logdir, "apply_warp2.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(
            f"applywarp -i {input.epi}  -o {output.epi2t1} -r {input.t1} --rel --interp=spline --premat={input.premat}"
        )
        shell(f"echo $(date) > {output.log}")


rule remove_bias_field:
    input:
        bias_field=f_BiasField,
        reg_epi=rules.apply_warp2.output.epi2t1,
    output:
        reg_epi_restore=join(outdir_T1w, "nodif2T1w_restore_initII.nii.gz"),
        log=join(logdir, "remove_bias_field.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(
            f"fslmaths {input.reg_epi} -div {input.bias_field} {output.reg_epi_restore}"
        )
        shell(f"echo $(date) > {output.log}")


rule bbregister:
    """
    Boundary-based registration of EPI to T1w.
    ===========================================
    The output of bbregister is a transformation matrix file (EPItoT1w.dat) that represents the spatial relationship between the EPI and T1-weighted images.
    """
    input:
        mov=rules.remove_bias_field.output.reg_epi_restore,
        init_reg=join(config["freesurferdir"], "{subid}", "mri/transforms/eye.dat"),
    output:
        nodif2t1=join(outdir_T1w, "nodif2T1w_bbreg.nii.gz"),
        reg=join(outdir_T1w, "EPItoT1w.dat"),
        log=join(logdir, "bbregister.done"),
    params:
        dof=dof,
        SUBJECTS_DIR=config["freesurferdir"],
        surf="white.deformed",
    group:
        "DiffusionToStructural"
    run:
        shell(
            f"""
            export SUBJECTS_DIR={params.SUBJECTS_DIR}
            bbregister --s {wildcards.subid} --mov {input.mov} --surf {params.surf} --init-reg {input.init_reg} --bold --reg {output.reg} --{params.dof} --o {output.nodif2t1}
            """
        )

        shell(f"echo $(date) > {output.log}")


rule tkregister2:
    """
    Fine-Tune registration
    ======================
    Uses a combination of intensity and surface-based registration to refine the EPI-to-T1w registration.
    The output of tkregister2 is a transformation matrix file (diff2str_fs.mat) that represents the spatial relationship between the EPI and T1-weighted images.
    """
    input:
        reg=rules.bbregister.output.reg,
        mov=rules.remove_bias_field.output.reg_epi_restore,
        targ=f_T1w,
    output:
        freesurfer_reg=join(outdir_T1w, "diff2str_fs.mat"),
        log=join(logdir, "tkregister2.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(
            f"tkregister2_cmdl --noedit --reg {input.reg} --mov {input.mov} --targ {input.targ} --fslregout {output.freesurfer_reg}"
        )
        shell(f"echo $(date) > {output.log}")


rule get_transformation_matrices:
    """
    The first command concatenates three transformation matrices: diff2str_fs.mat, "$regimg"2T1w_initII.mat, and diff2str.mat. The resulting transformation matrix is saved as diff2str.mat in the specified WorkingDirectory.
    The second command calculates the inverse of the diff2str.mat transformation matrix and saves it as str2diff.mat in the WorkingDirectory.
    """
    input:
        fsl_mat=rules.epi_reg_dof.output.mat,
        freesurfer_mat=rules.tkregister2.output.freesurfer_reg,
    output:
        diff2str=join(outdir_T1w, "diff2str.mat"),
        str2diff=join(outdir_T1w, "str2diff.mat"),
        log=join(logdir, "get_transformation_matrices.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(
            f"convert_xfm -omat {output.diff2str} -concat {input.freesurfer_mat} {input.fsl_mat}"
        )
        shell(f"convert_xfm -omat {output.str2diff} -inverse {output.diff2str}")
        shell(f"echo $(date) > {output.log}")


rule apply_warp3:
    input:
        epi=rules.mask_out_data.output.nodif,
        premat=rules.get_transformation_matrices.output.diff2str,
        t1=f_T1w,
    output:
        epi2t1=join(outdir_T1w, "nodif2T1w.nii.gz"),
        log=join(logdir, "apply_warp3.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(
            f"applywarp --rel --interp=spline -i {input.epi} -r {input.t1} -o {output.epi2t1} --premat={input.premat}"
        )
        shell(f"echo $(date) > {output.log}")


rule remove_bias:
    input:
        epi2t1=rules.apply_warp3.output.epi2t1,
        bias_field=f_BiasField,
    output:
        epi2t1_restore=join(outdir_T1w, "nodif2T1w_restore.nii.gz"),
        log=join(logdir, "remove_bias.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(
            f"fslmaths {input.epi2t1} -div {input.bias_field} {output.epi2t1_restore}"
        )
        shell(f"echo $(date) > {output.log}")


rule qc4:
    input:
        t1_restore=f_T1w_restore,
        epi2t1=rules.apply_warp3.output.epi2t1,
    output:
        nii=join(outdir_T1w, "QC_nodif.nii.gz"),
        fig=join(qcdir, "qc4.png"),
        log=join(logdir, "qc4.done"),
    group:
        "DiffusionToStructural"
    run:
        from nilearn import plotting
        import matplotlib.pyplot as plt
        import matplotlib

        shell(f"fslmaths {input.t1_restore} -mul {input.epi2t1} -sqrt {output.nii}")
        # Use a non-interactive backend
        matplotlib.use("agg")
        fig, axs = plt.subplots(nrows=2, figsize=(16, 9))
        plotting.plot_anat(output.nii, axes=axs[0])
        plotting.plot_stat_map(output.nii, bg_img=input.t1_restore, axes=axs[1])
        fig.savefig(output.fig)
        shell(f"echo $(date) > {output.log}")


rule diff_res_structural_space:
    """
    Generate a ${DiffRes} structural space for resampling the diffusion data into
    """
    input:
        data=rules.mask_out_data.output.data_masked,
        t1_restore=f_T1w_restore,
        epi2t1=rules.apply_warp3.output.epi2t1,
    output:
        t1_diffres=join(outdir_T1w, "T1w_acpc_dc_restore_diffRes.nii.gz"),
        diffres=join(outdir_T1w, "DiffRes.txt"),
        log=join(logdir, "diff_res_structural_space.done"),
    group:
        "DiffusionToStructural"
    run:
        import subprocess
        from nipype.interfaces.fsl import Info
        import subprocess

        cmd = f"fslval {input.data} pixdim1"
        result = subprocess.check_output(cmd, shell=True)
        DiffRes = eval(result.strip().decode("utf-8"))
        np.savetxt(output.diffres, [DiffRes], fmt="%0.2f")
        shell(
            f"flirt -interp spline -in {input.t1_restore} -ref {input.t1_restore} -applyisoxfm {DiffRes} -out {output.t1_diffres}"
        )
        shell(
            f"applywarp --rel --interp=spline -i {input.t1_restore} -r {output.t1_diffres} -o {output.t1_diffres}"
        )
        shell(f"echo $(date) > {output.log}")


rule diff_res_mask:
    input:
        fs_brain_mask=f_FreeSurferBrainMask,
        diffres=rules.diff_res_structural_space.output.diffres,
    output:
        t1_dwi_mask=join(outdir_T1w, "nodif_brain_mask.nii.gz"),
        log=join(logdir, "diff_res_mask.done"),
    group:
        "DiffusionToStructural"
    run:
        with open(input.diffres, "r") as f:
            DiffRes = f.read().strip()
        shell(
            f"flirt -interp nearestneighbour -in {input.fs_brain_mask} -ref {input.fs_brain_mask} -applyisoxfm {DiffRes} -out {output.t1_dwi_mask}"
        )
        shell(f"fslmaths {output.t1_dwi_mask} -kernel 3D -dilM {output.t1_dwi_mask}")
        shell(f"echo $(date) > {output.log}")


rule dilate_mask:
    input:
        mask=rules.diff_res_mask.output.t1_dwi_mask,
    output:
        dilated_mask=join(outdir_T1w, "nodif_brain_mask_temp.nii.gz"),
        log=join(logdir, "dilate_mask.done"),
    params:
        number_of_dilations=6,
    group:
        "DiffusionToStructural"
    run:
        shell(f"imcp {input.mask} {output.dilated_mask}")
        for i in range(params.number_of_dilations):
            shell(
                f"fslmaths {output.dilated_mask} -kernel 3D -dilM {output.dilated_mask}"
            )
        shell(f"echo $(date) > {output.log}")


rule rotate_bvecs_to_struct:
    input:
        bvals=rules.average_rotated_bvecs.output.av_bval,
        bvecs=rules.average_rotated_bvecs.output.av_bvec,
        mat=rules.get_transformation_matrices.output.diff2str,
    output:
        bvec_struct=join(outdir_T1w, "bvecs"),
        bvals_struct=join(outdir_T1w, "bvals"),
        log=join(logdir, "rotate_bvecs_to_struct.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(
            f"$HCPPIPEDIR/global/scripts/Rotate_bvecs.sh {input.bvecs} {input.mat} {output.bvec_struct}"
        )
        shell(f"cp {input.bvals} {output.bvals_struct}")
        shell(f"echo $(date) > {output.log}")


rule register_diff_to_T1w:
    # TODO: add option for Gradient Nonlinearity distortion correction
    input:
        data=rules.mask_out_data.output.data_masked,
        cnr_maps=rules.mask_out_data.output.cnr_maps_masked,
        fov_mask=rules.create_fov_mask.output.fov_mask,
        diffres=rules.diff_res_structural_space.output.diffres,
        t1_diffres=rules.diff_res_structural_space.output.t1_diffres,
        diff2str=rules.get_transformation_matrices.output.diff2str,
    output:
        data_dilated=join(datadir, "data_dilated.nii.gz"),
        cnr_maps_dilated=join(datadir, "cnr_maps_dilated.nii.gz"),
        data_t1=join(outdir_T1w, "data.nii.gz"),
        cnr_maps_t1=join(outdir_T1w, "cnr_maps.nii.gz"),
        fov_mask_t1=join(outdir_T1w, "fov_mask.nii.gz"),
        log=join(logdir, "register_diff_to_T1w.done"),
    group:
        "DiffusionToStructural"
    run:
        diffres = np.loadtxt(input.diffres, dtype=float)
        dilate_distance = diffres * 4
        # Dilation outside of the field of view to minimise the effect of the hard field of view edge on the interpolation
        shell(
            f"wb_command -volume-dilate {input.data} {dilate_distance} NEAREST {output.data_dilated}"
        )
        # Register diffusion data to T1w space without considering gradient nonlinearities
        shell(
            f"flirt -in {output.data_dilated} -ref {input.t1_diffres} -applyxfm -init {input.diff2str} -interp spline -out {output.data_t1}"
        )
        # Transform CNR maps
        shell(
            f"wb_command -volume-dilate {input.cnr_maps} {dilate_distance} NEAREST {output.cnr_maps_dilated}"
        )
        shell(
            f"flirt -in {output.cnr_maps_dilated} -ref {input.t1_diffres} -applyxfm -init {input.diff2str} -interp spline -out {output.cnr_maps_t1}"
        )
        # Transforms field of view mask to T1-weighted space
        shell(
            f"flirt -in {input.fov_mask} -ref {input.t1_diffres} -applyxfm -init {input.diff2str} -interp trilinear -out {output.fov_mask_t1}"
        )
        shell(f"echo $(date) > {output.log}")


rule adjust_fov_mask:
    """
    Only include voxels fully(!) within the field of view for every volume. This is done to avoid edge effects in the interpolation.
    """
    input:
        fov_mask=rules.register_diff_to_T1w.output.fov_mask_t1,
    output:
        fov_mask=join(outdir_T1w, "fov_mask_bin.nii.gz"),
        log=join(logdir, "adjust_fov_mask.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(f"fslmaths {input.fov_mask} -thr 0.999 -bin {output.fov_mask}")
        shell(f"echo $(date) > {output.log}")


rule fov_mask_data:
    input:
        fov_mask=rules.adjust_fov_mask.output.fov_mask,
        dilated_mask=rules.dilate_mask.output.dilated_mask,
        data=rules.register_diff_to_T1w.output.data_t1,
        cnr_maps=rules.register_diff_to_T1w.output.cnr_maps_t1,
    output:
        data=join(outdir_T1w, "data_masked.nii.gz"),
        cnr_maps=join(outdir_T1w, "QC", "cnr_maps.nii.gz"),
        log=join(logdir, "fov_mask_data.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(
            f"fslmaths {input.data} -mas {input.dilated_mask} -mas {input.fov_mask}  {output.data}"
        )
        shell(f"fslmaths {input.cnr_maps} -mas {input.fov_mask}  {output.cnr_maps}")
        shell(f"echo $(date) > {output.log}")


rule threshold_data:
    """
    Remove negative intensity values (from eddy) from final data
    """
    input:
        data=rules.fov_mask_data.output.data,
    output:
        data=join(outdir_T1w, "data_masked_thresh.nii.gz"),
        log=join(logdir, "threshold_data.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(f"fslmaths {input.data} -thr 0 {output.data}")
        shell(f"echo $(date) > {output.log}")


rule new_nodif_brain_mask:
    input:
        mask=rules.diff_res_mask.output.t1_dwi_mask,
        data=rules.threshold_data.output.data,
    output:
        tmean=join(outdir_T1w, "_temp_tmean.nii.gz"),
        mask=join(outdir_T1w, "nodif_brain_mask_new.nii.gz"),
        log=join(logdir, "new_nodif_brain_mask.done"),
    group:
        "DiffusionToStructural"
    run:
        shell(f"fslmaths {input.data} -Tmean {output.tmean}")
        shell(f"fslmaths {input.mask} -mas {output.tmean} {output.mask}")
        shell(f"echo $(date) > {output.log}")


rule calculate_coverage:
    input:
        nodif_brain_mask_old=rules.diff_res_mask.output.t1_dwi_mask,
        nodif_brain_mask=rules.new_nodif_brain_mask.output.mask,
    output:
        stats=join(outdir_T1w, "nodif_brain_mask.stats.txt"),
        log=join(logdir, "calculate_coverage.done"),
    group:
        "DiffusionToStructural"
    run:
        import subprocess

        cmd1 = f"fslstats {input.nodif_brain_mask_old} -V | awk '{{print $1}}'"
        NvoxBrainMask = (
            subprocess.check_output(cmd1, shell=True).decode("utf-8").strip()
        )
        cmd2 = f"fslstats {input.nodif_brain_mask} -V | awk '{{print $1}}'"
        NvoxFinalMask = (
            subprocess.check_output(cmd2, shell=True).decode("utf-8").strip()
        )
        cmd3 = f"echo 'scale=4; 100 * {NvoxFinalMask} / {NvoxBrainMask}' | bc -l"
        PctCoverage = subprocess.check_output(cmd3, shell=True).decode("utf-8").strip()
        shell(f"echo 'PctCoverage, NvoxFinalMask, NvoxBrainMask' >> {output.stats}")
        shell(
            f"echo '{PctCoverage}, {NvoxFinalMask}, {NvoxBrainMask}' >> {output.stats}"
        )
        shell(f"echo $(date) > {output.log}")


################################
# MRTrix3                      #
################################

MRTrix_out = join(outdir_T1w, "MRTrix")
MRTrix_log = join(MRTrix_out, "log")


rule create_mif:
    """Create .mif for further processing with MRTrix."""
    input:
        data=rules.threshold_data.output.data,
        bvals=rules.rotate_bvecs_to_struct.output.bvals_struct,
        bvecs=rules.rotate_bvecs_to_struct.output.bvec_struct,
    output:
        dwi=join(MRTrix_out, "DWI.mif"),
        log=join(logdir, "create_mif.done"),
    log:
        stdout=MRTrix_log,
        sterr=MRTrix_log,
    group:
        "preprocessing"
    threads: 1
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    run:
        shell(
            f"mrconvert {input.data} {output.dwi} --fslgrad {input.bvecs} {input.bvals} -datatype float32 -stride 0,0,0,1 -nthreads {threads} -force -info"
        )
        shell(f"echo $(date) > {output.log}")


rule tissue_segmentation:
    """Run FSL-based tissue segmentation algorithm"""
    input:
        t1_brain=join(dir_T1w, "T1w_acpc_dc_restore_brain.nii.gz"),
        t2_brain=join(dir_T1w, "T2w_acpc_dc_restore_brain.nii.gz"),
    output:
        seg=join(MRTrix_out, "5TT.mif"),
        log=join(logdir, "tissue_segmentation.done"),
    group:
        "preprocessing"
    threads: 16
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    run:
        shell(
            f"5ttgen fsl {input.t1_brain} {output.seg} -t2 {input.t2_brain} -premasked -nthreads {threads}"
        )
        shell(f"echo $(date) > {output.log}")


rule bias_correction:
    """ANTS-based bias correction."""
    input:
        dwi=rules.create_mif.output.dwi,
    output:
        dwi_bias=join(MRTrix_out, "DWI_bias_ants.mif"),
        dwi_bias_field=join(MRTrix_out, "DWI_bias_ants_field.mif"),
        log=join(logdir, "bias_correction.done"),
    group:
        "preprocessing"
    threads: 1
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    container:
        config["containers"]["mrtrix"]
    shell:
        r"""
        dwibiascorrect ants {input.dwi} {output.dwi_bias} \
            -bias {output.dwi_bias_field} \
            -nthreads {threads} \
            -force -info

        echo $(date) > {output.log}
        """


rule create_mean_b0:
    """Mean image"""
    input:
        dwi_bias=rules.bias_correction.output.dwi_bias,
    output:
        mean_b0=join(MRTrix_out, "meanb0.mif"),
        log=join(logdir, "create_mean_b0.done"),
    group:
        "preprocessing"
    threads: 1
    resources:
        mem_mb=8 * 1024,
        runtime=5,  # in minutes
    run:
        shell(
            f"dwiextract {input.dwi_bias} - -bzero | mrmath - mean {output.mean_b0} -axis 3 -force -info"
        )
        shell(f"echo $(date) > {output.log}")


rule extract_response_functions:
    """Extract response function. Uses -stride 0,0,0,1"""
    input:
        dwi_bias=rules.bias_correction.output.dwi_bias,
    output:
        response_wm=join(MRTrix_out, "response_wm.txt"),
        response_gm=join(MRTrix_out, "response_gm.txt"),
        response_csf=join(MRTrix_out, "response_csf.txt"),
        voxels=join(MRTrix_out, "RF_voxels.mif"),
        log=join(logdir, "extract_response_functions.done"),
    group:
        "preprocessing"
    threads: 1
    # container:
    #     config["containers"]["mrtrix"]
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    run:
        shell(
            f"dwi2response dhollander {input.dwi_bias} {output.response_wm} {output.response_gm} {output.response_csf} -voxels {output.voxels} -nthreads {threads} -force -info"
        )
        shell(f"echo $(date) > {output.log}")


rule generate_mask:
    """Create DWI brain mask,"""
    input:
        dwi_bias=join(MRTrix_out, "DWI_bias_ants.mif"),
    output:
        dwi_mask=join(MRTrix_out, "DWI_mask.mif"),
    group:
        "preprocessing"
    threads: 1
    resources:
        mem_mb=8 * 1024,
        runtime=10,  # in minutes
    run:
        shell(
            f"dwi2mask {input.dwi_bias} {output.dwi_mask} -nthreads {threads} -force -info"
        )
