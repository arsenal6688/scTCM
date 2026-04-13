from sctcm import recode4anno

# 配置路径
ADATA_PATH = "/scTCM/sctcm/tests/result/sc_marker4anno-marker_annotated.h5ad"
ANNOTATION_TEMPLATE = "/scTCM/docs/templates/sc_recode4anno-input_template.txt"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# 执行分析流程
if __name__ == "__main__":
    adata = recode4anno.run_annotation_pipeline(
        adata_path=ADATA_PATH,
        template_path=ANNOTATION_TEMPLATE,
        output_dir=OUTPUT_DIR,
        cluster_key="leiden_res1_0"
    )

print("✅ recode4anno 测试流程执行完成！")
