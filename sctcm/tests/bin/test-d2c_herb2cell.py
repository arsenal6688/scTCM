from sctcm.d2c import run_herb2cell

# ====================== 路径 ======================
SC_TEMPLATE = "/scTCM/docs/templates/sc_io-input_template.txt"
ANNDATA_PATH = "/scTCM/sctcm/tests/result/sc_recode4anno-celltype_manual_annotated.h5ad"
HERB_PKL = "/scTCM/sctcm/tests/result/tcm_herb2target-herb2target.pkl"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

if __name__ == "__main__":
    run_herb2cell(
        sc_template=SC_TEMPLATE,
        adata_path=ANNDATA_PATH,
        herb_pkl=HERB_PKL,
        output_dir=OUTPUT_DIR,
        min_pct=0.25,
        logfc_threshold=2.0,
        use_filtered_only=False   
    )

print("✅ herb2cell 运行完成！")
