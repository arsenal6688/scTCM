# -*- coding: utf-8 -*-
import os
import pickle

def run_form2nonnested(
    input_pkl: str,
    output_dir: str
):
    """
    将三层嵌套的 tcm_form2targetall.pkl 转换为一层扁平结构
    原结构：{方剂: {中药: {化合物: [基因...]}}}
    新结构：{方剂: [所有基因...去重]}
    """
    os.makedirs(output_dir, exist_ok=True)
    output_pkl = os.path.join(output_dir, "tcm_form2nonnested-non_nested.pkl")

    # 加载原始三层嵌套数据
    with open(input_pkl, "rb") as f:
        data = pickle.load(f)

    # 扁平化：方剂 → 所有基因（去重）
    nonnested_dict = {}
    for form_name, herb_dict in data.items():
        all_genes = []
        # 遍历方剂下的每味中药
        for herb_name, compound_dict in herb_dict.items():
            # 遍历中药下的每个化合物
            for gene_list in compound_dict.values():
                all_genes.extend(gene_list)
        # 去重 + 排序
        nonnested_dict[form_name] = sorted(list(set(all_genes)))

    # 保存
    with open(output_pkl, "wb") as f:
        pickle.dump(nonnested_dict, f)

    print(f"✅ 方剂非嵌套转换完成")
    print(f"📂 输出文件：{output_pkl}")
    print(f"🧾 方剂总数：{len(nonnested_dict)}")

    return nonnested_dict
