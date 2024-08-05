import os
import nibabel as nib
import argparse
from hcpy.utils import extract_b0_indices, get_pedir, get_readout_time

def extract_b0_volumes(input_nii, input_bval, output_b0, b0maxbval, b0dist):
    # Get b0 indices
    indices, log = extract_b0_indices(input_bval, b0maxbval, b0dist=b0dist)

    indcount = 1  # Initialize index counter

    # Prepare FSL command to extract b0 volumes
    for idx in indices:
        output_file = os.path.splitext(output_b0)[0] + "_{}.nii.gz".format(idx)
        os.system("fslroi {} {} {} 1".format(input_nii, output_file, idx))
        indcount += 1

    PEdir = get_pedir(input_nii)
    ro_time = get_readout_time(input_nii, PEdir)

    # Merge the extracted b0 volumes
    output_b0_dir = os.path.splitext(output_b0)[0]
    shell_cmd = "fslmerge -t {} {}_*.nii.gz".format(output_b0, output_b0_dir)
    os.system(shell_cmd)

    # Remove temporary files
    output_b0_base = os.path.splitext(output_b0)[0]
    shell_cmd = "rm {}_*.nii.gz".format(output_b0_base)
    os.system(shell_cmd)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract b0 volumes from DWI data.')
    parser.add_argument('--input_nii', type=str, required=True, help='Input NIfTI file')
    parser.add_argument('--input_bval', type=str, required=True, help='Input bval file')
    parser.add_argument('--output_b0', type=str, required=True, help='Output b0 file')
    parser.add_argument('--b0maxbval', type=float, required=True, help='Maximum b-value for b0 volumes')
    parser.add_argument('--b0dist', type=float, required=True, help='Distance between b0 volumes')

    args = parser.parse_args()

    extract_b0_volumes(args.input_nii, args.input_bval, args.output_b0, args.b0maxbval, args.b0dist)