from sctcm.tcm import run_formula2target

TEMPLATE_PATH = "/scTCM/docs/templates/tcm_formula2target-template.txt"
FORMULA_BROWSE_PATH = "/scTCM/sctcm/data/BATMAN-TCM/4-Formulas_information_and_their_herbal_components/formula_browse.txt"
HERB_BROWSE_PATH = "/scTCM/sctcm/data/BATMAN-TCM/3-Herbs_information_and_their_ingredient_compounds/herb_browse.txt"
HERB_INFO_PATH = "/scTCM/sctcm/data/HERB_ingredient_info_v2.txt"
KNOWN_TARGET_PATH = "/scTCM/sctcm/data/BATMAN-TCM/1-TCM_ingredient_compounds_information_and_their_protein_targets/known_browse_by_ingredients.txt.gz"
OUTPUT_PKL = "/scTCM/sctcm/tests/result/tcm_formula2target-formula2target.pkl"

# ====================== 执行 ======================
if __name__ == "__main__":
    formula2target_dict = run_formula2target(
        template_path=TEMPLATE_PATH,
        formula_browse_path=FORMULA_BROWSE_PATH,
        herb_browse_path=HERB_BROWSE_PATH,
        herb_info_path=HERB_INFO_PATH,
        known_target_path=KNOWN_TARGET_PATH,
        output_pkl=OUTPUT_PKL,
        drug_likeness_min=0.1,
        ob_score_min=10.0,
        use_filtered_only=False
    )

print("✅ formula2target 测试流程执行完成！")
print(f"✅ 共处理 {len(formula2target_dict)} 个方剂")
