import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")))

from sctcm.tcm.form2nonnested import run_form2nonnested

# ===================== 路径配置 =====================
INPUT_PKL = "/scTCM/sctcm/tests/result/tcm_form2targetall.pkl"
OUTPUT_DIR = "/scTCM/sctcm/tests/result"

# ====================== 执行 ======================
if __name__ == "__main__":
    print("=== 开始转换：form2targetall → 非嵌套结构 ===")
    
    form_target_nested = run_form2nonnested(
        input_pkl=INPUT_PKL,
        output_dir=OUTPUT_DIR
    )

print("\n🎉 全部完成！")
print(f"✅ 最终非嵌套结构：{len(form_target_nested)} 个方剂")
