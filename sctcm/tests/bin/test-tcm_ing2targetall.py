import sys
sys.path.insert(0, "/scTCM")

from sctcm.tcm.ing2targetall import run_ing2targetall

HERB_PATH = "/scTCM/sctcm/data/HERB_ingredient_info_v2.txt"
KNOWN_PATH = "/scTCM/sctcm/data/BATMAN-TCM/1-TCM_ingredient_compounds_information_and_their_protein_targets/known_browse_by_ingredients.txt.gz"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

if __name__ == "__main__":
    print("="*60)
    result = run_ing2targetall(
        herb_info_path=HERB_PATH,
        known_target_path=KNOWN_PATH,
        output_dir=OUTPUT_DIR,
        ob_min=0.0,
        drug_likeness_min=0.0
    )
    print(f"✅ 最终成功：{len(result)} 个成分")
    print("="*60)
