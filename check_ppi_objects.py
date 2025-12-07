"""
简化的PPI距离检查脚本

在PyMOL中运行后，手动检查几个distance对象的距离
"""

from pymol import cmd

print("\n" + "="*60)
print("PPI Distance对象列表")
print("="*60)

# 获取所有ppi_开头的distance对象
all_objects = cmd.get_names("objects")
ppi_objects = [obj for obj in all_objects if obj.startswith("ppi_")]

if not ppi_objects:
    print("\n⚠️ 没有找到PPI distance对象")
else:
    print(f"\n找到 {len(ppi_objects)} 个PPI distance对象:")
    
    # 按类型分组
    by_type = {}
    for obj in ppi_objects:
        # 提取类型 (ppi_hbond_0 -> hbond)
        parts = obj.split("_")
        if len(parts) >= 2:
            obj_type = parts[1]
            if obj_type not in by_type:
                by_type[obj_type] = []
            by_type[obj_type].append(obj)
    
    print("\n按类型分组:")
    for obj_type, objs in sorted(by_type.items()):
        print(f"\n  {obj_type}: {len(objs)} 个")
        # 只显示前3个
        for obj in objs[:3]:
            print(f"    - {obj}")
        if len(objs) > 3:
            print(f"    ... 还有 {len(objs)-3} 个")

print("\n" + "="*60)
print("颜色说明:")
print("="*60)
print("  hbond (氢键)       -> 蓝色 (glue_blue)")
print("  saltbridge (盐桥)  -> 红色 (glue_red)")
print("  hydrophobic (疏水) -> 绿色 (glue_green)")
print("  pipi (π-π堆积)     -> 黄色 (yellow)")
print("  cationpi (阳离子-π)-> 紫色 (purple)")
print("="*60)

print("\n💡 提示:")
print("   1. 在PyMOL中点击任意虚线，查看底部显示的距离")
print("   2. 或者运行: get_distance ppi_hbond_0")
print("   3. 正常的相互作用距离应该 <5Å")
print("="*60 + "\n")
