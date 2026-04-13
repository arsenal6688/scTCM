# sctcm/__init__.py
__version__ = "0.1.0"  # 版本号，用户可通过sctcm.__version__查看

# 导出核心子模块（用户可直接import sctcm.sc / sctcm.tcm）
from . import sc
from .sc import dimred
from .sc import cluster
from .sc import marker4anno
from .sc import AUCell4anno
from .sc import recode4anno

from . import tcm
from . import d4c
from . import pl
from . import utils


# 定义__all__：明确*导入时的暴露范围（规范）
__all__ = ["sc", "tcm", "d4c", "pl", "utils", "__version__"]
