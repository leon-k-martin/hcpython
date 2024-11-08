import numpy as np
from nipype.interfaces.fsl import Merge
import nibabel as nib
from hcpy import utils
import argparse
import json

def merge_data(imgs, b0_imgs, bvals, bvecs, b0maxbval, b0dist, output):
    np.set_printoptions(precision=6)

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

    for img, b0_img, bval, bvec in zip(imgs.split(), b0_imgs.split(), bvals.split(), bvecs.split()):
        series_index = np.append(series_index, np.repeat(ind, nib.load(img).shape[3]))

        acqparams, b0_volumes, log = utils.get_b0_acqparms(img, bval, b0maxbval=b0maxbval, b0dist=b0dist)
        bval = np.loadtxt(bval)

        pedir = utils.get_pedir(img)
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

    with open(output['b0_log'], "w") as f:
        f.write("\n".join(pos_log + neg_log))

    with open(output['topup_acqparams'], "w") as f:
        f.write("\n".join(pos_acqparams + neg_acqparams))

    with open(output['eddy_acqparams'], "w") as f:
        f.write("\n".join(pos_acqparams + neg_acqparams))

    np.savetxt(output['series_index'], series_index, fmt="%i")

    pos_b0s_str = " ".join(pos_b0s)
    neg_b0s_str = " ".join(neg_b0s)
    posneg_imgs = " ".join(pos_imgs + neg_imgs)

    merge_pos_b0 = Merge(in_files=pos_b0s, dimension='t', merged_file=output['b0_pos'])
    merge_pos_b0.run()

    merge_neg_b0 = Merge(in_files=neg_b0s, dimension='t', merged_file=output['b0_neg'])
    merge_neg_b0.run()

    merge_posneg_b0 = Merge(in_files=[output['b0_pos'], output['b0_neg']], dimension='t', merged_file=output['b0_posneg'])
    merge_posneg_b0.run()

    merge_posneg_imgs = Merge(in_files=pos_imgs + neg_imgs, dimension='t', merged_file=output['img_posneg'])
    merge_posneg_imgs.run()

    np.savetxt(output['bval_pos'], pos_bvals.reshape(-1)[None], fmt="%i", delimiter=" ")
    np.savetxt(output['bval_neg'], neg_bvals.reshape(-1)[None], fmt="%i", delimiter=" ")
    np.savetxt(output['bvec_pos'], pos_bvecs, fmt="%.6f", delimiter=" ")
    np.savetxt(output['bvec_neg'], neg_bvecs, fmt="%.6f", delimiter=" ")

    bvals = np.append(pos_bvals, neg_bvals, axis=0)
    bvecs = np.append(pos_bvecs, neg_bvecs, axis=1)

    index = utils.get_eddy_index(bvals, b0maxbval=b0maxbval, b0dist=b0dist)
    np.savetxt(output['index'], index, fmt="%i")

    np.savetxt(output['bvals_posneg'], bvals.reshape(-1)[None], fmt="%i", delimiter=" ")
    np.savetxt(output['bvecs_posneg'], bvecs, fmt="%.6f", delimiter=" ")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge neuroimaging data.")
    parser.add_argument("--imgs", type=str, required=True, help="Input image file paths separated by spaces")
    parser.add_argument("--b0_imgs", type=str, required=True, help="Input b0 image file paths separated by spaces")
    parser.add_argument("--bvals", type=str, required=True, help="Input b-values file paths separated by spaces")
    parser.add_argument("--bvecs", type=str, required=True, help="Input b-vectors file paths separated by spaces")
    parser.add_argument("--b0maxbval", type=float, required=True, help="Maximum b-value for b0 volumes")
    parser.add_argument("--b0dist", type=float, required=True, help="Distance between b0 volumes")
    parser.add_argument("--output", type=str, required=True, help="Output file paths as JSON string")

    args = parser.parse_args()

    output_paths = json.loads(args.output)

    merge_data(args.imgs, args.b0_imgs, args.bvals, args.bvecs, args.b0maxbval, args.b0dist, output_paths)