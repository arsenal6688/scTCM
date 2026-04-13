import os
import scanpy as sc
import omicverse as ov
import matplotlib.pyplot as plt
import pandas as pd
from typing import Optional, List, Union, Dict
from anndata import AnnData
import warnings
warnings.filterwarnings("ignore")

# 定义默认参数
DEFAULT_SUBTYPE_PARAMS = {
    "neighbors": {
        "n_neighbors": 15,
        "n_pcs": 30,
        "use_rep": "X_harmony",
        "key_added": None
    },
    "leiden": {
        "resolutions": 1.0,
        "key_prefix": "subtype_leiden_res",
        "neighbors_key": None
    },
    "mde": {
        "use_rep": "X_harmony",
        "mde_key": "X_mde_subtype"
    },
    "plot": {
        "basis": "X_mde_subtype",
        "dpi": 300,
        "frameon": "small",
        "palette": None
    },
    "preprocess": {
        "re_scale": False,
        "layer": "scaled"
    }
}

def read_cell_type_template(template_path: str) -> List[str]:
    if not os.path.exists(template_path):
        raise FileNotFoundError(f"输入模板文件不存在：{template_path}")
    cell_types = []
    with open(template_path, "r", encoding="utf-8") as f:
        for line in f.readlines():
            ct = line.strip()
            if ct:
                cell_types.append(ct)
    if not cell_types:
        raise ValueError("输入模板文件未检测到有效细胞类型")
    cell_types = list(set(cell_types))
    print(f"从模板中提取到待分析的主细胞类型：{cell_types}")
    return cell_types

def extract_major_cell_type(adata: AnnData, cell_type: str, cell_type_key: str = "cell_type") -> AnnData:
    if cell_type_key not in adata.obs.columns:
        common_cell_type_keys = ["celltype_manual", "cell_type", "CellType", "celltype"]
        matched_keys = [k for k in common_cell_type_keys if k in adata.obs.columns]
        if matched_keys:
            cell_type_key = matched_keys[0]
            print(f"自动匹配到细胞类型列：{cell_type_key}")
        else:
            raise KeyError(f"未找到细胞类型列：{cell_type_key}，可用列：{adata.obs.columns.tolist()}")

    if cell_type not in adata.obs[cell_type_key].unique():
        raise ValueError(f"细胞类型 {cell_type} 不存在，可用类型：{adata.obs[cell_type_key].unique()}")

    adata_sub = adata[adata.obs[cell_type_key] == cell_type, :].copy()
    print(f"提取 {cell_type} 完成：{adata.n_obs} → {adata_sub.n_obs} 个细胞")
    adata_sub.obs = adata_sub.obs.reset_index(drop=True)
    
    # 打印提取细胞类型的基础统计信息
    print(f"📊 {cell_type} 基础统计：")
    print(f"   - 总细胞数：{adata_sub.n_obs}")
    print(f"   - 总基因数：{adata_sub.n_vars}")
    if 'n_counts' in adata_sub.obs.columns:
        print(f"   - 平均每个细胞的UMI数：{adata_sub.obs['n_counts'].mean():.2f}")
    if 'n_genes' in adata_sub.obs.columns:
        print(f"   - 平均每个细胞的基因数：{adata_sub.obs['n_genes'].mean():.2f}")
    
    return adata_sub

def print_subtype_statistics(adata_sub: AnnData, cluster_key: str, resolution: float):
    """打印亚型的详细统计信息"""
    subtype_counts = adata_sub.obs[cluster_key].value_counts().sort_index()
    total_cells = len(adata_sub)
    
    print(f"\n📈 分辨率 {resolution} 亚型统计详情：")
    print(f"   总计：{len(subtype_counts)} 个亚型，{total_cells} 个细胞")
    print("   ┌─────────────┬───────────┬────────────┐")
    print("   │  亚型ID     │  细胞数量  │  占比(%)   │")
    print("   ├─────────────┼───────────┼────────────┤")
    
    for subtype, count in subtype_counts.items():
        percentage = (count / total_cells) * 100
        print(f"   │  {subtype:<9} │  {count:>7}  │  {percentage:>8.2f} │")
    
    print("   └─────────────┴───────────┴────────────┘")
    print(f"   最大亚型：{subtype_counts.idxmax()} ({subtype_counts.max()} 个细胞, {subtype_counts.max()/total_cells*100:.2f}%)")
    print(f"   最小亚型：{subtype_counts.idxmin()} ({subtype_counts.min()} 个细胞, {subtype_counts.min()/total_cells*100:.2f}%)")
    print(f"   平均每个亚型细胞数：{subtype_counts.mean():.2f}")
    print(f"   亚型细胞数标准差：{subtype_counts.std():.2f}")
    
    # 保存统计信息到adata的uns中
    if 'subtype_statistics' not in adata_sub.uns:
        adata_sub.uns['subtype_statistics'] = {}
    adata_sub.uns['subtype_statistics'][cluster_key] = {
        'counts': subtype_counts.to_dict(),
        'total_cells': total_cells,
        'mean_cells_per_subtype': subtype_counts.mean(),
        'std_cells_per_subtype': subtype_counts.std(),
        'max_subtype': subtype_counts.idxmax(),
        'max_subtype_count': subtype_counts.max(),
        'min_subtype': subtype_counts.idxmin(),
        'min_subtype_count': subtype_counts.min()
    }

def run_subtype_analysis(adata_sub: AnnData, params: Dict = None) -> AnnData:
    final_params = DEFAULT_SUBTYPE_PARAMS.copy()
    if params is not None:
        for k, v in params.items():
            if isinstance(v, dict) and k in final_params:
                final_params[k].update(v)

    if final_params["preprocess"]["re_scale"]:
        ov.pp.scale(adata_sub, layer=final_params["preprocess"]["layer"])

    print("\n===== 构建亚型邻域图 =====")
    sc.pp.neighbors(
        adata_sub,
        n_neighbors=final_params["neighbors"]["n_neighbors"],
        n_pcs=final_params["neighbors"]["n_pcs"],
        use_rep=final_params["neighbors"]["use_rep"]
    )
    print(f"   邻域图构建完成：使用 {final_params['neighbors']['use_rep']} 作为表达矩阵，n_neighbors={final_params['neighbors']['n_neighbors']}, n_pcs={final_params['neighbors']['n_pcs']}")

    print("\n===== 执行亚型Leiden聚类 =====")
    resolutions = final_params["leiden"]["resolutions"]
    if not isinstance(resolutions, list):
        resolutions = [resolutions]

    for res in resolutions:
        res_key = f"{final_params['leiden']['key_prefix']}{str(res).replace('.', '_')}"
        sc.tl.leiden(
            adata_sub,
            resolution=res,
            key_added=res_key,
            random_state=42
        )
        print(f"\n分辨率 {res}：识别 {len(adata_sub.obs[res_key].unique())} 个亚型")
        
        # 打印详细的亚型统计信息
        print_subtype_statistics(adata_sub, res_key, res)
        
        adata_sub.obs["subtype_default"] = adata_sub.obs[res_key]

    print("\n===== 执行亚型MDE降维 =====")

    adata_sub.obsm[final_params["mde"]["mde_key"]] = ov.utils.mde(
        adata_sub.obsm[final_params["mde"]["use_rep"]]
    )

    print(f"MDE降维完成：生成 {final_params['mde']['mde_key']} 维度 {adata_sub.obsm[final_params['mde']['mde_key']].shape[1]}")
    
    # 打印整体分析统计摘要
    print("\n📋 亚型分析整体统计摘要：")
    print(f"   - 输入细胞总数：{adata_sub.n_obs}")
    print(f"   - 聚类使用分辨率：{resolutions}")
    print(f"   - 最终识别亚型数：{len(adata_sub.obs['subtype_default'].unique())}")
    print(f"   - MDE降维维度：{adata_sub.obsm[final_params['mde']['mde_key']].shape}")
    
    return adata_sub

def plot_subtype_cluster(adata_sub: AnnData, cell_type: str, save_dir: str, params=None):
    final_params = DEFAULT_SUBTYPE_PARAMS.copy()
    if params:
        for k, v in params.items():
            if isinstance(v, dict) and k in final_params:
                final_params[k].update(v)

    cluster_keys = [k for k in adata_sub.obs.columns if k.startswith(final_params["leiden"]["key_prefix"])]
    for ck in cluster_keys:
        res = ck.replace(final_params["leiden"]["key_prefix"], "").replace("_", ".")
        save_path = os.path.join(save_dir, f"sc_subtypeProc-{cell_type.lower()}_subtype_res{res.replace('.','_')}.png")
        ov.utils.embedding(adata_sub, basis=final_params["plot"]["basis"], color=ck, show=False)
        plt.savefig(save_path, dpi=final_params["plot"]["dpi"], bbox_inches="tight")
        plt.close()
        print(f"图片已保存：{save_path}")

def save_subtype_h5ad(adata_sub, cell_type, save_dir):
    h5ad_path = os.path.join(save_dir, f"sc_subtypeProc-{cell_type.lower()}.h5ad")
    adata_sub.write(h5ad_path, compression="gzip")
    print(f"h5ad已保存：{h5ad_path}")
    
    # 保存统计信息到CSV文件
    stats_path = os.path.join(save_dir, f"sc_subtypeProc-{cell_type.lower()}_statistics.csv")
    all_stats = []
    if 'subtype_statistics' in adata_sub.uns:
        for res_key, stats in adata_sub.uns['subtype_statistics'].items():
            res = res_key.replace(DEFAULT_SUBTYPE_PARAMS['leiden']['key_prefix'], '').replace('_', '.')
            for subtype, count in stats['counts'].items():
                all_stats.append({
                    'resolution': res,
                    'subtype_id': subtype,
                    'cell_count': count,
                    'percentage': (count / stats['total_cells']) * 100,
                    'total_cells': stats['total_cells'],
                    'mean_cells_per_subtype': stats['mean_cells_per_subtype'],
                    'std_cells_per_subtype': stats['std_cells_per_subtype']
                })
    if all_stats:
        pd.DataFrame(all_stats).to_csv(stats_path, index=False, encoding='utf-8')
        print(f"统计信息已保存到CSV：{stats_path}")
    
    return h5ad_path

def run_subtype_proc_pipeline(
    adata_path,
    template_path,
    save_root,
    cell_type_key="cell_type",
    custom_params=None
):
    os.makedirs(save_root, exist_ok=True)
    cell_types = read_cell_type_template(template_path)
    adata = ov.read(adata_path)
    
    # 打印整体数据集信息
    print(f"\n📁 输入数据集信息：")
    print(f"   - 总细胞数：{adata.n_obs}")
    print(f"   - 总基因数：{adata.n_vars}")
    print(f"   - 细胞类型列：{cell_type_key}")
    print(f"   - 唯一细胞类型数：{len(adata.obs[cell_type_key].unique())}")
    
    res = {}
    summary_stats = []

    for ct in cell_types:
        print(f"\n{'='*50}\n开始分析 {ct} 亚型\n{'='*50}")
        sub = extract_major_cell_type(adata, ct, cell_type_key)
        sub = run_subtype_analysis(sub, custom_params)
        plot_subtype_cluster(sub, ct, save_root, custom_params)
        h5ad = save_subtype_h5ad(sub, ct, save_root)
        res[ct] = h5ad
        
        # 收集汇总统计信息
        n_subtypes = len(sub.obs['subtype_default'].unique())
        summary_stats.append({
            'cell_type': ct,
            'cell_count': sub.n_obs,
            'n_subtypes': n_subtypes,
            'h5ad_path': h5ad
        })

    # 打印最终汇总统计
    print(f"\n{'='*60}")
    print("📊 所有细胞类型亚型分析汇总")
    print("="*60)
    summary_df = pd.DataFrame(summary_stats)
    print(summary_df.to_string(index=False))
    print(f"\n总计分析细胞类型数：{len(summary_df)}")
    print(f"总计分析细胞数：{summary_df['cell_count'].sum()}")
    print(f"总计识别亚型数：{summary_df['n_subtypes'].sum()}")
    
    # 保存汇总统计
    summary_path = os.path.join(save_root, "sc_subtypeProc-subtype_analysis_summary.csv")
    summary_df.to_csv(summary_path, index=False, encoding='utf-8')
    print(f"\n汇总统计已保存：{summary_path}")
    
    print("\n✅ 全部分析完成！")
    return res
