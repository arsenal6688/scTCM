from sctcm.d4c.TCMscore2cell import TCMscore2cell

# 路径
INPUT_TXT = "/scTCM/docs/templates/d4c_TCMscore2cell-input_template.txt"
SC_TXT = "/scTCM/docs/templates/sc_io-input_template.txt"
ANNO_H5AD = "/scTCM/sctcm/tests/result/sc_recode4anno-celltype_manual_annotated.h5ad"
FORM_CONFIG = "/scTCM/docs/templates/d4c_TCMscore2cell-formula_herb_config_template.txt"
OUT = "/scTCM/sctcm/tests/result"

if __name__ == "__main__":
    TCMscore2cell(
        input_txt=INPUT_TXT,
        sc_meta_txt=SC_TXT,
        anno_h5ad=ANNO_H5AD,
        formula_config_txt=FORM_CONFIG,
        output_dir=OUT
    )
