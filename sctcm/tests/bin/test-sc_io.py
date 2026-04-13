import os
import sys

from sctcm.sc.io import read_and_merge_samples, write_h5ad

INPUT_TEMPLATE = "/scTCM/docs/templates/sc_io-input_template.txt"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"
OUTPUT_H5AD = os.path.join(OUTPUT_DIR, "sc_io-merged_samples.h5ad")

adata_merged = read_and_merge_samples(
        sample_table=INPUT_TEMPLATE,
        sep="\t"
)

write_h5ad(adata_merged, OUTPUT_H5AD)

