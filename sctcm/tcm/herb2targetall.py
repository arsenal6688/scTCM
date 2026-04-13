#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd
import pickle
import re
import warnings
warnings.filterwarnings("ignore")

__all__ = ["run_herb2targetall"]

def parse_herb_ingredients(herb_browse_path):
    print("读取 herb_browse.txt ...")
    herb_df = pd.read_csv(herb_browse_path, sep="\t", encoding="utf-8")
    herb_df = herb_df.dropna(subset=["Ingredients"]).copy()

    def extract_cids(ing_str):
        ids = []
        for ing in ing_str.split("|"):
            m = re.search(r"\((\d+)\)", ing)
            if m:
                ids.append(m.group(1))
        return ids

    herb_df["cids"] = herb_df["Ingredients"].apply(extract_cids)
    herb_df = herb_df[herb_df["cids"].apply(len) > 0].copy()
    print(f"有效中药数：{len(herb_df)}")
    return herb_df

def load_ing2target(pkl_path):
    with open(pkl_path, "rb") as f:
        return pickle.load(f)

def build_herb2target_nested(herb_df, ing2target):
    herb2target = {}

    for _, row in herb_df.iterrows():
        pinyin = str(row["Pinyin.Name"]).strip()
        cid_list = row["cids"]

        component_dict = {}
        for cid in cid_list:
            targets = ing2target.get(cid, [])
            if len(targets) > 0:
                component_dict[cid] = targets

        if len(component_dict) > 0:
            herb2target[pinyin] = component_dict

    print(f"最终有靶点的中药数：{len(herb2target)}")
    return herb2target

def run_herb2targetall(herb_browse_path, ing_target_pkl, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    herb_df = parse_herb_ingredients(herb_browse_path)
    ing2target = load_ing2target(ing_target_pkl)
    herb2target = build_herb2target_nested(herb_df, ing2target)

    out_pkl = os.path.join(output_dir, "tcm_herb2targetall.pkl")
    with open(out_pkl, "wb") as f:
        pickle.dump(herb2target, f)

    print(f"\n✅ 已保存：{out_pkl}")
    return herb2target
