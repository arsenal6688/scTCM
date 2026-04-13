from sctcm.tcm import run_ingredient2target

# ====================== 配置路径 ======================
# 输入模板文件
TEMPLATE_PATH = "/scTCM/docs/templates/tcm_ingredient2target-template.txt"
# HERB成分信息文件
HERB_INFO_PATH = "/scTCM/sctcm/data/HERB_ingredient_info_v2.txt"
# 已知靶点文件
KNOWN_TARGET_PATH = "/scTCM/sctcm/data/BATMAN-TCM/1-TCM_ingredient_compounds_information_and_their_protein_targets/known_browse_by_ingredients.txt.gz"
# 输出文件目录
OUTPUT_PKL = "/scTCM/sctcm/tests/result/tcm_ingredient2target-ingredient2target.pkl"

# ====================== 执行分析流程 ======================
if __name__ == "__main__":
    cid2target_dict = run_ingredient2target(
        template_path=TEMPLATE_PATH,
        herb_info_path=HERB_INFO_PATH,
        known_target_path=KNOWN_TARGET_PATH,
        output_pkl=OUTPUT_PKL,
        drug_likeness_min=0.1,
        ob_score_min=10.0,
        use_filtered_only=False
    )

print("✅ ingredient2target 测试流程执行完成！")
print(f"✅ 共生成 {len(cid2target_dict)} 个成分的靶点映射")
