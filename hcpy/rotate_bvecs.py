# rotate_bvecs.py
import sys
import numpy as np
import os
import subprocess
from nipype.interfaces import fsl
import numpy as np


def main(input_bvecs, matrix_file, output_bvecs):
    # Check for file existence
    if not os.path.exists(input_bvecs):
        print(f"Source bvecs {input_bvecs} does not exist!")
        sys.exit(1)
    if not os.path.exists(matrix_file):
        print(f"Matrix file {matrix_file} does not exist!")
        sys.exit(1)

    # Extract rotation matrix from the affine
    rotation_matrix = extract_rotation_matrix(matrix_file)

    # Load bvecs
    bvecs = np.loadtxt(input_bvecs)

    # Ensure bvecs is 3xn
    if bvecs.shape[0] != 3:
        bvecs = bvecs.T

    # Apply rotation
    rotated_bvecs = np.dot(rotation_matrix, bvecs)

    # Save rotated bvecs
    np.savetxt(output_bvecs, rotated_bvecs, fmt="%10.6f")


def extract_rotation_matrix(matrix_file):
    """
    Extracts the rotation matrix from a FLIRT affine transformation file using Nipype's interface to FSL's avscale.

    Parameters:
    - matrix_file: Path to the FLIRT affine matrix file.

    Returns:
    - A 3x3 numpy array containing the rotation matrix.
    """
    # Initialize AvScale with the path to the matrix file
    avscale = fsl.AvScale()
    avscale.inputs.mat_file = matrix_file

    # Execute avscale
    try:
        result = avscale.run()
    except Exception as e:
        print("Error running AvScale with Nipype:", e)
        return None

    # Extract the rotation matrix from the result
    # The rotation matrix is part of the 'rotations' key in the avscale output
    # Convert the list of lists to a NumPy array
    rotation_matrix = np.array(result.outputs.rotation)

    return rotation_matrix


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: <original bvecs> <affine matrix> <rotated (output) bvecs>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
