# ======================== 测试脚本 ========================
# 路径参数
SC_META_TXT = "/scTCM/docs/templates/sc_io-input_template.txt"
ANNDATA_PATH = "/scTCM/sctcm/tests/result/d4c_TCMscoreScan-annotated.h5ad"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# ======================== 执行 ========================
from sctcm.pl.TCMscoreScan_case_high import TCMscoreScan_case_high

if __name__ == "__main__":
    TCMscoreScan_case_high(
        adata_path=ANNDATA_PATH,
        output_dir=OUTPUT_DIR,
        group_col="Group",
        case_label="CASE",
        ctrl_label="CTRL",
        celltype_col="celltype_manual",
        pval_cutoff=0.05,
        fc_cutoff=1.5
    )
