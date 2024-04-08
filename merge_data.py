def merge_data(input, params, output):
    import subprocess
    import numpy as np
    import nibabel as nib
    from hcpy.utils import get_pedir, get_eddy_index, get_b0_acqparms, get_readout_time

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
        series_index = np.append(series_index, np.repeat(ind, nib.load(img).shape[3]))

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

    subprocess.run(f"fslmerge -t {output.b0_pos} {' '.join(pos_b0s)}")
    subprocess.run(f"fslmerge -t {output.b0_neg} {' '.join(neg_b0s)}")
    subprocess.run(f"fslmerge -t {output.b0_posneg} {output.b0_pos} {output.b0_neg}")

    subprocess.run(f"fslmerge -t {output.img_posneg} {' '.join(pos_imgs + neg_imgs)}")

    np.savetxt(output.bval_pos, pos_bvals.reshape(-1)[None], fmt="%i", delimiter=" ")
    np.savetxt(output.bval_neg, neg_bvals.reshape(-1)[None], fmt="%i", delimiter=" ")
    np.savetxt(output.bvec_pos, pos_bvecs, fmt="%.6f", delimiter=" ")
    np.savetxt(output.bvec_neg, neg_bvecs, fmt="%.6f", delimiter=" ")

    bvals = np.append(pos_bvals, neg_bvals, axis=0)
    bvecs = np.append(pos_bvecs, neg_bvecs, axis=1)

    # Create index file
    index = get_eddy_index(bvals, b0maxbval=params.b0maxbval, b0dist=params.b0dist)
    np.savetxt(output.index, index, fmt="%i")

    # Save bvals and bvecs
    np.savetxt(output.bvals_posneg, bvals.reshape(-1)[None], fmt="%i", delimiter=" ")
    np.savetxt(output.bvecs_posneg, bvecs, fmt="%.6f", delimiter=" ")
