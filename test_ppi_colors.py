"""
PPI颜色测试脚本 - 验证自定义颜色是否正确注册

在PyMOL中运行:
run /Users/vesper/Desktop/git/GlueTK/GlueTK/test_ppi_colors.py
"""

from pymol import cmd

print("\n" + "="*60)
print("测试PPI自定义颜色注册")
print("="*60)

# 注册自定义颜色
glue_red = [215/255.0, 92/255.0, 93/255.0]
glue_blue = [64/255.0, 124/255.0, 174/255.0]
glue_green = [142/255.0, 186/255.0, 141/255.0]

cmd.set_color("glue_red", glue_red)
cmd.set_color("glue_blue", glue_blue)
cmd.set_color("glue_green", glue_green)

print("\n✅ 自定义颜色已注册:")
print(f"  glue_red   = RGB({int(glue_red[0]*255)}, {int(glue_red[1]*255)}, {int(glue_red[2]*255)})")
print(f"  glue_blue  = RGB({int(glue_blue[0]*255)}, {int(glue_blue[1]*255)}, {int(glue_blue[2]*255)})")
print(f"  glue_green = RGB({int(glue_green[0]*255)}, {int(glue_green[1]*255)}, {int(glue_green[2]*255)})")

# 验证颜色是否可用
try:
    # 创建测试对象
    cmd.pseudoatom("test_red", pos=[0, 0, 0])
    cmd.pseudoatom("test_blue", pos=[5, 0, 0])
    cmd.pseudoatom("test_green", pos=[10, 0, 0])
    
    # 应用颜色
    cmd.color("glue_red", "test_red")
    cmd.color("glue_blue", "test_blue")
    cmd.color("glue_green", "test_green")
    
    # 创建测试distance对象
    cmd.distance("test_dist_red", "test_red", "test_blue")
    cmd.distance("test_dist_green", "test_blue", "test_green")
    
    # 设置虚线颜色
    cmd.color("glue_red", "test_dist_red")
    cmd.set("dash_color", "glue_red", "test_dist_red")
    cmd.set("dash_width", 3.0, "test_dist_red")
    
    cmd.color("glue_green", "test_dist_green")
    cmd.set("dash_color", "glue_green", "test_dist_green")
    cmd.set("dash_width", 3.0, "test_dist_green")
    
    print("\n✅ 测试对象已创建:")
    print("  - test_red (红色球)")
    print("  - test_blue (蓝色球)")
    print("  - test_green (绿色球)")
    print("  - test_dist_red (红色虚线)")
    print("  - test_dist_green (绿色虚线)")
    print("\n👀 请检查PyMOL视图中的颜色是否正确")
    print("   如果虚线仍然是黄色，说明PyMOL版本可能不支持自定义dash_color")
    
except Exception as e:
    print(f"\n❌ 测试失败: {e}")

print("="*60 + "\n")
