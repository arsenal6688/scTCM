from sctcm.tcm import run_herb2target

TEMPLATE_PATH = "/scTCM/docs/templates/tcm_herb2target-template.txt"
HERB_BROWSE_PATH = "/scTCM/sctcm/data/BATMAN-TCM/3-Herbs_information_and_their_ingredient_compounds/herb_browse.txt"
HERB_INFO_PATH = "/scTCM/sctcm/data/HERB_ingredient_info_v2.txt"
KNOWN_TARGET_PATH = "/scTCM/sctcm/data/BATMAN-TCM/1-TCM_ingredient_compounds_information_and_their_protein_targets/known_browse_by_ingredients.txt.gz"
OUTPUT_PKL = "/scTCM/sctcm/tests/result/tcm_herb2target-herb2target.pkl"

# ====================== 执行 ======================
if __name__ == "__main__":
    herb2target_dict = run_herb2target(
        template_path=TEMPLATE_PATH,
        herb_browse_path=HERB_BROWSE_PATH,
        herb_info_path=HERB_INFO_PATH,
        known_target_path=KNOWN_TARGET_PATH,
        output_pkl=OUTPUT_PKL,
        drug_likeness_min=0.1,
        ob_score_min=10.0,
        use_filtered_only=False
    )

print("✅ herb2target 测试流程执行完成！")
print(f"✅ 共处理 {len(herb2target_dict)} 味中药")

