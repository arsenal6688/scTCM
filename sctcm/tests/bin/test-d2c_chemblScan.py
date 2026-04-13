# -*- coding: utf-8 -*-
from sctcm.d2c import chemblScan

# ====================== 配置路径 ======================
# 1. 筛选完成的 Figure2k 轻量 h5ad 文件（含B细胞、内皮、动脉/静脉/心内膜）
ANNDATA_PATH = "/script4scTCMmethodologyPaper/01.ingredientCompa/case1/fig2k_light.h5ad"

# 2. ChEMBL36 药物-靶点字典 PKL
CHEMBL_PKL = "/scTCM/sctcm/tests/result/chembl_36_nested_drug_dict.pkl"

# 3. 输出结果文件夹
OUTPUT_DIR = "/script4scTCMmethodologyPaper/01.ingredientCompa/case1"

# 4. 细胞类型分组列名
GROUP_BY = "cell_type_plot"

# ====================== 执行分析 ======================
if __name__ == "__main__":
    chemblScan.run_chembl_drug2cell(
        h5ad_path=ANNDATA_PATH,
        drug_dict_pkl=CHEMBL_PKL,
        output_dir=OUTPUT_DIR,
        groupby=GROUP_BY,
        n_processes=8,
        top_n=30,
        plot_top=15
    )
