from sctcm.d2c import run_formula2cell

# ====================== 路径 ======================
SC_TEMPLATE = "/scTCM/docs/templates/sc_io-input_template.txt"
ANNDATA_PATH = "/scTCM/sctcm/tests/result/sc_recode4anno-celltype_manual_annotated.h5ad"
FORMULA_PKL = "/scTCM/sctcm/tests/result/tcm_formula2target-formula2target.pkl"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

if __name__ == "__main__":
    run_formula2cell(
        sc_template=SC_TEMPLATE,
        adata_path=ANNDATA_PATH,
        formula_pkl=FORMULA_PKL,
        output_dir=OUTPUT_DIR,
        min_pct=0.25,
        logfc_threshold=2.0,
        use_filtered_only=False   
    )

print("✅ formula2cell 运行完成！")
