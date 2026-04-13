from sctcm.d2c import formulaScan

# ====================== 配置路径 ======================
SC_TEMPLATE = "/scTCM/docs/templates/sc_io-input_template.txt"
ANNDATA_PATH = "/scTCM/sctcm/tests/result/sc_recode4anno-celltype_manual_annotated.h5ad"
FORMULA_PKL = "/scTCM/sctcm/tests/result/tcm_form2targetall.pkl"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# ====================== 执行 ======================
if __name__ == "__main__":
    formulaScan(
        sc_template=SC_TEMPLATE,
        adata_path=ANNDATA_PATH,
        formula_pkl=FORMULA_PKL,
        output_dir=OUTPUT_DIR,
        min_pct=0.05,
        logfc_threshold=2.0,
        use_filtered_only=True
    )

print("✅ formulaScan 全复方批量drug2cell分析运行完成！")
