from sctcm.d2c import herbScan

# ====================== 配置路径======================
# 1. 单细胞模板文件（获取Group分组信息）
SC_TEMPLATE = "/scTCM/docs/templates/sc_io-input_template.txt"

# 2. 手动注释完成的h5ad文件
ANNDATA_PATH = "/scTCM/sctcm/tests/result/sc_recode4anno-celltype_manual_annotated.h5ad"

# 3. 全草药-靶点字典 PKL 文件
HERB_PKL = "/scTCM/sctcm/tests/result/tcm_herb2targetall.pkl"

# 4. 输出结果文件夹
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# ====================== 执行分析 ======================
if __name__ == "__main__":
    herbScan(
        sc_template=SC_TEMPLATE,
        adata_path=ANNDATA_PATH,
        herb_pkl=HERB_PKL,
        output_dir=OUTPUT_DIR,
        min_pct=0.1,
        logfc_threshold=2.0,
        use_filtered_only=True
    )

print("✅ herbScan 全草药批量drug2cell分析运行完成！")
