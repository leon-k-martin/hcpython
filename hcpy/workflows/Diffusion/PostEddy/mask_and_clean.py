from nipype import Workflow, Node
from nipype.interfaces.fsl import ExtractROI, BET
from nipype.interfaces.fsl.maths import MathsCommand, ApplyMask, Threshold
from nipype.interfaces.utility import IdentityInterface
import os
from nipype.interfaces.io import DataSink


def setup_workflow(input, params, output):
    # Define the workflow
    wf = Workflow(name="mask_and_clean_workflow")
    wf.base_dir = params.logdir

    # DataSink setup
    datasink = Node(DataSink(), name="datasink")
    datasink.inputs.base_directory = params.datadir
    datasink.inputs.container = ""

    # Define an input node
    inputspec = Node(
        IdentityInterface(fields=["data", "fov_mask", "cnr_maps"]), name="inputspec"
    )
    inputspec.inputs.data = input.data
    inputspec.inputs.fov_mask = input.fov_mask
    inputspec.inputs.cnr_maps = input.cnr_maps

    # Mask out data outside the field of view for 'data' and 'cnr_maps'
    mask_data = Node(ApplyMask(), name="mask_data", out_file=input.data)
    mask_cnr_maps = Node(ApplyMask(), name="mask_cnr_maps", out_file=input.cnr_maps)

    # Remove negative intensity values from 'data'
    remove_negative = Node(Threshold(thresh=0), name="remove_negative")

    # Extract the first volume (assuming this is the operation intended by fslroi with 0 1)
    extract_nodif = Node(
        ExtractROI(t_min=0, t_size=1), name="extract_nodif", out_file=output.nodif
    )

    # Brain extraction on the first volume
    bet_nodif = Node(
        BET(frac=0.1, mask=True), name="bet_nodif", out_file=output.nodif_brain
    )
    wf.connect(
        [
            (
                inputspec,
                mask_data,
                [("data", "in_file"), ("fov_mask", "mask_file")],
            ),
            (
                inputspec,
                mask_cnr_maps,
                [("cnr_maps", "in_file"), ("fov_mask", "mask_file")],
            ),
            (mask_data, remove_negative, [("out_file", "in_file")]),
            (remove_negative, extract_nodif, [("out_file", "in_file")]),
            (extract_nodif, bet_nodif, [("roi_file", "in_file")]),
        ]
    )

    # Output connections to DataSink
    wf.connect(remove_negative, "out_file", datasink, "@data")
    wf.connect(mask_cnr_maps, "out_file", datasink, "@cnr_maps")
    wf.connect(extract_nodif, "roi_file", datasink, "@nodif")
    wf.connect(bet_nodif, "out_file", datasink, "@nodif_brain")
    wf.connect(bet_nodif, "mask_file", datasink, "@nodif_brain_mask")
    return wf
