# ====================== 配置路径  ======================
# 1. 扫描输入模板（成分/单味药/复方 配置）
SCAN_TEMPLATE = "/scTCM/docs/templates/d4c_TCMscoreScan-input_template.txt"

# 2. 单细胞模板文件（获取Group分组信息）
SC_TEMPLATE = "/scTCM/docs/templates/sc_io-input_template.txt"

# 3. 手动注释完成的h5ad文件
ANNDATA_PATH = "/scTCM/sctcm/tests/result/sc_recode4anno-celltype_manual_annotated.h5ad"

# 4. 输出结果文件夹
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# ====================== 执行分析 ======================
from sctcm.d4c import TCMscoreScan

if __name__ == "__main__":
    TCMscoreScan(
        scan_input_txt=SCAN_TEMPLATE,
        sc_meta_txt=SC_TEMPLATE,
        anno_h5ad=ANNDATA_PATH,
        output_dir=OUTPUT_DIR
    )
