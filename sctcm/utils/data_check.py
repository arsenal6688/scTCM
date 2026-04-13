# sctcm/utils/data_check.py
import anndata as ad
#from collections import dict

def check_anndata(adata: ad.AnnData, required_obsm: list = None, required_obs: list = None) -> bool:
    """
    校验AnnData对象是否包含必要的字段
    :param adata: AnnData对象
    :param required_obsm: 必需的obsm字段（如["X_umap"]）
    :param required_obs: 必需的obs字段（如["celltype"]）
    :return: 校验通过返回True，否则抛异常
    """
    if not isinstance(adata, ad.AnnData):
        raise TypeError("输入必须是AnnData对象")
    # 校验obsm
    if required_obsm is not None:
        for key in required_obsm:
            if key not in adata.obsm:
                raise ValueError(f"AnnData缺少obsm字段：{key}")
    # 校验obs
    if required_obs is not None:
        for key in required_obs:
            if key not in adata.obs.columns:
                raise ValueError(f"AnnData缺少obs字段：{key}")
    return True

def check_tcm_dict(tcm_dict: dict) -> bool:
    """
    校验中医药字典格式（key为字符串，value为列表）
    """
    if not isinstance(tcm_dict, dict):
        raise TypeError("输入必须是字典")
    for k, v in tcm_dict.items():
        if not isinstance(k, str):
            raise ValueError(f"字典key必须是字符串：{k}")
        if not isinstance(v, list):
            raise ValueError(f"字典value必须是列表：{k}")
    return True
