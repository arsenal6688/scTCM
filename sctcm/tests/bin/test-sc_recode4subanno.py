#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from sctcm.sc import recode4anno

# ===================== 【路径配置】 =====================
# 1. 巨噬细胞亚型h5ad文件（subtypeProc输出）
ADATA_PATH = "/scTCM/sctcm/tests/result/sc_subtypeProc-macrophage.h5ad"
# 2. 亚型注释模板（0~11 对应具体亚型名称）
ANNOTATION_TEMPLATE = "/scTCM/docs/templates/sc_recode4anno-subtype_input_template.txt"
# 3. 结果输出目录
OUTPUT_DIR = "/scTCM/sctcm/tests/result/Macrophage"  

# ===================== 【核心参数】 =====================
CLUSTER_KEY = "subtype_leiden_res1_0"  
NEW_ANNO_KEY = "subtype_manual"        

# ===================== 【一键执行亚型重新注释流程】 =====================
if __name__ == "__main__":
    print(f"🚀 开始巨噬细胞亚型手动注释...")
    print(f"📋 配置信息：")
    print(f"  - 输入文件：{ADATA_PATH}")
    print(f"  - 注释模板：{ANNOTATION_TEMPLATE}")
    print(f"  - 聚类列名：{CLUSTER_KEY}")
    print(f"  - 输出目录：{OUTPUT_DIR}")
    
    # 核心调用：复用recode4anno，仅修改参数适配亚型
    adata_annotated = recode4anno.run_annotation_pipeline(
        adata_path=ADATA_PATH,
        template_path=ANNOTATION_TEMPLATE,
        output_dir=OUTPUT_DIR,
        cluster_key=CLUSTER_KEY,       # 指定亚型聚类列
        new_anno_key=NEW_ANNO_KEY,     # 新增亚型注释列
        run_umap=False,                
        run_mde=True                   # 绘制MDE图展示注释结果（与subtypeProc可视化一致）
    )
    
    # 打印注释结果统计
    print(f"\n✅ 亚型注释完成！")
    print(f"📊 注释结果统计：")
    print(adata_annotated.obs[NEW_ANNO_KEY].value_counts())
    print(f"\n📁 输出文件清单：")
    print(f"  - 注释后h5ad：{OUTPUT_DIR}/sc_recode4anno-celltype_manual_annotated.h5ad")
    print(f"  - MDE可视化图：{OUTPUT_DIR}/sc_recode4anno-manual_annotation_mde.png")
    print(f"  - 备份文件：{OUTPUT_DIR}/sc_recode4anno-backup_uns.pkl")
