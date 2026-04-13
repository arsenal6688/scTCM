from sctcm.sc import subtypeProc

# 配置路径
ADATA_PATH = "/scTCM/sctcm/tests/result/sc_recode4anno-celltype_manual_annotated.h5ad"
INPUT_TEMPLATE = "/scTCM/docs/templates/sc_subtypeProc-input_template.txt"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# 执行分析流程
if __name__ == "__main__":
    result_h5ad_dict = subtypeProc.run_subtype_proc_pipeline(
        adata_path=ADATA_PATH,
        template_path=INPUT_TEMPLATE,
        save_root=OUTPUT_DIR,
        cell_type_key="celltype_manual",
        custom_params={
            "neighbors": {"n_neighbors": 15, "n_pcs": 30},
            "leiden": {"resolutions": 1.0}
        }
    )

    print("✅ 亚型分析测试流程执行完成！")
