# -*- coding: utf-8 -*-
import os
import pickle

def run_herb2nonnested(
    input_pkl: str,
    output_dir: str
):
    """
    将嵌套结构的 tcm_herb2targetall.pkl 转换为非嵌套扁平结构
    原结构：{草药: {化合物ID: [基因...]}}
    新结构：{草药: [所有基因...去重]}

    Parameters
    ----------
    input_pkl : str
        原始嵌套 herb2targetall.pkl 路径
    output_dir : str
        输出目录

    Returns
    -------
    dict
        非嵌套草药-靶点字典 {草药名: [基因1, 基因2...]}
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    output_pkl = os.path.join(output_dir, "tcm_herb2nonnested-non_nested.pkl")


    # 加载原始嵌套数据
    with open(input_pkl, "rb") as f:
        data = pickle.load(f)

    # 转换为扁平结构
    nonnested_dict = {}
    for herb_name, compound_target_dict in data.items():
        all_genes = []
        # 收集所有基因
        for gene_list in compound_target_dict.values():
            all_genes.extend(gene_list)
        # 去重并保存
        nonnested_dict[herb_name] = sorted(list(set(all_genes)))

    # 保存
    with open(output_pkl, "wb") as f:
        pickle.dump(nonnested_dict, f)

    print(f"✅ 转换完成：非嵌套草药靶点字典已保存")
    print(f"📂 输出路径：{output_pkl}")
    print(f"🌿 草药总数：{len(nonnested_dict)}")

    return nonnested_dict



__all__ = ["run_herb2nonnested"]
