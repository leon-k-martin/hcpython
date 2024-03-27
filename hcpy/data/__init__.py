from os.path import dirname, realpath, join
import pandas as pd

DATA_DIR = dirname(realpath(__file__))

fs_lut = join(DATA_DIR, "FreeSurferColorLUT.txt")

mmp_annot_lh = join(DATA_DIR, "lh.HCP-MMP1.annot")
mmp_annot_rh = join(DATA_DIR, "rh.HCP-MMP1.annot")

mmp_mrtrix_original = join(DATA_DIR, "hcpmmp1_original.txt")
mmp_mrtrix_ordered = join(DATA_DIR, "hcpmmp1_ordered.txt")

mmp_atlas_info = pd.read_csv(
    mmp_mrtrix_ordered,
    comment="#",
    sep="\s+",
    header=None,
    names=["ID", "Labelname", "R", "G", "B", "A"],
)

dk_mrtrix_ordered = join(DATA_DIR, "fs_default.txt")
dk_atlas_info = pd.read_csv(
    dk_mrtrix_ordered,
    comment="#",
    sep="\s+",
    header=None,
    names=["ID", "Acronym", "Labelname", "R", "G", "B", "A"],
)
