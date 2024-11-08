import json
import nibabel as nib
import os
import argparse
from hcpy.utils import get_pedir


def create_series_file(
    input_imgs, output_pos_corr, output_neg_corr, output_pos, output_neg
):
    pos_vols = []
    neg_vols = []

    with open(output_pos_corr, "w") as pos_corr_file, open(
        output_neg_corr, "w"
    ) as neg_corr_file, open(output_pos, "w") as pos_file, open(
        output_neg, "w"
    ) as neg_file:

        for img in input_imgs.split(" "):
            volumes = nib.load(img).shape[3]
            pedir = get_pedir(img)
            if "-" in pedir:
                neg_vols.append(volumes)
            else:
                pos_vols.append(volumes)

        pos_vols.sort()
        neg_vols.sort()

        for vol_pos, vol_neg in zip(pos_vols, neg_vols):
            min_vol = min(vol_pos, vol_neg)
            pos_file.write("{} {}\n".format(min_vol, vol_pos))
            neg_file.write("{} {}\n".format(min_vol, vol_neg))
            if vol_pos > 0:
                pos_corr_file.write("{}\n".format(min_vol))
            if vol_neg > 0:
                neg_corr_file.write("{}\n".format(min_vol))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create series file from neuroimaging data."
    )
    parser.add_argument(
        "--input_imgs",
        type=str,
        required=True,
        help="Input image file paths separated by spaces",
    )
    parser.add_argument(
        "--output_pos_corr",
        type=str,
        required=True,
        help="Output positive correlation file",
    )
    parser.add_argument(
        "--output_neg_corr",
        type=str,
        required=True,
        help="Output negative correlation file",
    )
    parser.add_argument(
        "--output_pos", type=str, required=True, help="Output positive file"
    )
    parser.add_argument(
        "--output_neg", type=str, required=True, help="Output negative file"
    )

    args = parser.parse_args()

    create_series_file(
        args.input_imgs,
        args.output_pos_corr,
        args.output_neg_corr,
        args.output_pos,
        args.output_neg,
    )
