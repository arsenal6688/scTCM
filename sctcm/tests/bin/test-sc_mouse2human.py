import os
import sys

from sctcm.sc.io import read_and_merge_samples
from sctcm.sc.mouse2human import process_mouse_sc_data

INPUT_TEMPLATE = "/scTCM/docs/templates/sc_mouse2human-input_template.txt"
HOMOLOG_FILE = "/scTCM/sctcm/data/mouse_human_homologs_dec2021.csv"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"
OUTPUT_H5AD = os.path.join(OUTPUT_DIR, "sc_mouse2human-converted.h5ad")

# 执行小鼠转人基因流程
process_mouse_sc_data(
    sample_table=INPUT_TEMPLATE,
    homolog_file=HOMOLOG_FILE,
    output_h5ad=OUTPUT_H5AD,
    sep="\t",
    min_cells=1
)
