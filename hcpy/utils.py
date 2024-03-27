import os
import json
import numpy as np
from os.path import join, basename, dirname
import hcpy
import subprocess


def get_pedir(img):
    with open(img.replace(".nii.gz", ".json"), "r") as j:
        metadata = json.load(j)
    pedir = metadata["PhaseEncodingDirection"]
    return pedir


def get_readout_time(f_nii, PEdir):
    """
    Calculates the total readout time for an MRI scan based on its JSON metadata and NIfTI file.

    Parameters:
    f_nii (str): The file path to the NIfTI image. If a JSON file path is provided instead,
                 it will be adjusted to the corresponding NIfTI file path.
    PEdir (str): The phase encoding direction, which can be specified with characters
                 indicating the direction (e.g., 'i', 'j', 'L', 'P').

    Returns:
    float: The total readout time in seconds.

    The function first identifies the correct JSON metadata file associated with the
    NIfTI image. It then reads this JSON file to extract the EffectiveEchoSpacing, which
    is essential for calculating the readout time. The function also determines the
    dimension of the image in the phase encoding direction (dimP) by invoking the `fslval`
    command-line utility from FSL. Using the EffectiveEchoSpacing and the dimension size,
    it calculates the readout time considering all phase encoding steps except the first one
    (hence dimP - 1).

    The readout time is finally converted to seconds before being returned.

    Note: This function requires the subprocess module to call `fslval` for obtaining
    dimension sizes and assumes the presence of a corresponding JSON metadata file for the
    NIfTI image.
    """
    if not f_nii.endswith(".json"):
        f_json = f_nii.replace(".nii.gz", ".json")
    else:
        f_json = f_nii
        f_nii = f_nii.replace(".json", ".nii.gz")
    with open(f_json) as json_file:
        metadata = json.load(json_file)
        # Extract the echo spacing from the metadata
    echospacing = metadata.get("EffectiveEchoSpacing", 0)
    # Determine the dimension based on the phase encoding direction
    if "i" in PEdir or "L" in PEdir:
        dimP = int(subprocess.check_output([f"fslval", f_nii, "dim1"]))
    elif "j" in PEdir or "P" in PEdir:
        dimP = int(subprocess.check_output([f"fslval", f_nii, "dim2"]))
    dimPminus1 = dimP - 1
    # Calculate the total readout time
    ro_time_ms = float(echospacing) * dimPminus1
    # ro_time_secs = ro_time_ms / 1000
    # Write the total readout time to the output file
    return ro_time_ms


def get_acqparams(f_nii):
    PEdir = get_pedir(f_nii)
    ro_time = get_readout_time(f_nii, PEdir)
    if "-" in PEdir:
        value = -1
    else:
        value = 1

    if "j" in PEdir:
        acqparams = f"0 {value} 0 {ro_time:6f}"
    elif "i" in PEdir:
        acqparams = f"{value} 0 0 {ro_time:6f}"
    return acqparams


# def extract_b0_indices(bval_file, b0maxbval):
#     with open(bval_file, "r") as f:
#         bvals = [int(bval) for bval in f.readline().strip().split()]

#     return [i for i, bval in enumerate(bvals) if bval < b0maxbval]


def extract_b0_indices(bval_file, b0maxbval=50, b0dist=45):
    if isinstance(bval_file, str):
        bvals = np.loadtxt(bval_file, dtype=int)
    elif isinstance(bval_file, np.ndarray):
        bvals = bval_file
    else:
        raise ValueError("bval_file must be a string or numpy array")
    # with open(bval_file, "r") as f:
    #     bvals = [int(bval) for bval in f.readline().strip().split()]
    log = ""
    b0_indices = []
    # Initialize to a value that ensures the first b0 can be selected if it meets criteria:
    last_b0_index = -(b0dist + 10)

    for i, bval in enumerate(bvals):
        if bval < b0maxbval and (i - last_b0_index) > b0dist:
            b0_indices.append(i)
            last_b0_index = i  # Update the index of the last b0 volume added
            log += f"Extracting Volume {i} as a b=0. Measured b={bval}\n"
    return b0_indices, log


def get_eddy_index(bval, b0maxbval=50, b0dist=45):
    if isinstance(bval, str):
        bval = np.loadtxt(bval, dtype=int)

    b0_volumes, log = extract_b0_indices(bval, b0maxbval=b0maxbval, b0dist=b0dist)
    b0_ind = 1
    index = np.zeros(len(bval))
    for b0 in b0_volumes:
        index[b0:] = b0_ind
        b0_ind += 1
    return index


def get_b0_acqparms(f_nii, bval, b0maxbval=50, b0dist=45):
    b0_indices, log = extract_b0_indices(bval, b0maxbval, b0dist)
    acqparams = []
    for i in b0_indices:
        acqparams.append(get_acqparams(f_nii))
    return acqparams, b0_indices, log


def get_sift_threshold(sift_weights, n_tck=200000):
    if isinstance(sift_weights, str):
        sift_weights = np.loadtxt(sift_weights)
    # Calculate the index for the kth largest element
    k = len(sift_weights) - n_tck  # 99,800,000
    # Use numpy.partition to find the kth largest weight
    # The kth largest will be at index k in the partitioned array
    threshold = np.partition(sift_weights, k)[k]
    return threshold
