import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")))

from sctcm.ChEMBL.parseChEMBLfiltered import run_parse_chembl36

# ===================== 修改阈值 =====================
MY_THRESHOLDS = {
    'none': 6,        # 1 μM
    'NHR': 7,         # 100 nM
    'GPCR': 7,        # 100 nM
    'Ion Channel': 5, # 10 μM
    'Kinase': 7.53,   # 30 nM
}
# ==============================================================

INPUT_CHEMBL_DB = "/scTCM/sctcm/data/chembl_36_sqlite/chembl_36.db"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

print("=== Starting ChEMBL36 parsing for drug2cell ===")
run_parse_chembl36(INPUT_CHEMBL_DB, OUTPUT_DIR, thresholds_dict=MY_THRESHOLDS)
print("=== All done successfully ===")
