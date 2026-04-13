import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")))

# 导入转换函数
from sctcm.tcm.herb2nonnested import run_herb2nonnested

# ===================== 路径配置 =====================
INPUT_PKL = "/scTCM/sctcm/tests/result/tcm_herb2targetall.pkl"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# ====================== 执行 ======================
if __name__ == "__main__":
    print("=== 开始转换：herb2targetall → 非嵌套结构 ===")
    
    herb_target_nested = run_herb2nonnested(
        input_pkl=INPUT_PKL,
        output_dir=OUTPUT_DIR
    )

print("\n🎉 全部完成！")
print(f"✅ 最终非嵌套结构：{len(herb_target_nested)} 味草药")
