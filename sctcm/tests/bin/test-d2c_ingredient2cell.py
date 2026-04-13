from sctcm.d2c import run_ingredient2cell

# ====================== 配置路径 ======================
# 1. 单细胞模板（获取Group名称）
SC_TEMPLATE = "/scTCM/docs/templates/sc_io-input_template.txt"

# 2. 已注释好的h5ad路径
ANNDATA_PATH = "/scTCM/sctcm/tests/result/sc_recode4anno-celltype_manual_annotated.h5ad"

# 3. 成分-靶点PKL
INGREDIENT_PKL = "/scTCM/sctcm/tests/result/tcm_ingredient2target-ingredient2target.pkl"

# 4. 输出文件夹
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# ====================== 执行 ======================
if __name__ == "__main__":
    run_ingredient2cell(
        sc_template=SC_TEMPLATE,
        adata_path=ANNDATA_PATH,
        ingredient_pkl=INGREDIENT_PKL,
        output_dir=OUTPUT_DIR,
        min_pct=0.25,
        logfc_threshold=2.0,
        use_filtered_only=False  
    )

print("✅ ingredient2cell 运行完成！")
