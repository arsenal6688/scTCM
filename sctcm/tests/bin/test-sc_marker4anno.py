from sctcm import marker4anno

# 配置路径
ADATA_PATH = "/scTCM/sctcm/tests/result/sc_cluster.h5ad"
MARKER_FILE = "/scTCM/docs/templates/sc_marker4anno-marker_dict.txt"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# 执行分析流程
if __name__ == "__main__":
    adata = marker4anno.marker_analysis_pipeline(
        adata_path=ADATA_PATH,
        marker_file_path=MARKER_FILE,
        output_dir=OUTPUT_DIR,
        groupby="leiden_res1_0",
        deg_method="wilcoxon",
        top_n_degs=10,
        dpi=300
    )
print("✅ 测试流程执行完成！")
