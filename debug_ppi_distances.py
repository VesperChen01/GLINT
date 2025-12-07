"""
PPI距离调试脚本 - 检查distance对象的实际距离值

在PyMOL中运行:
run /Users/vesper/Desktop/git/GlueTK/GlueTK/debug_ppi_distances.py
"""

from pymol import cmd

print("\n" + "="*60)
print("检查PPI Distance对象")
print("="*60)

# 获取所有ppi_开头的distance对象
all_objects = cmd.get_names("objects")
ppi_objects = [obj for obj in all_objects if obj.startswith("ppi_")]

if not ppi_objects:
    print("\n⚠️ 没有找到PPI distance对象")
    print("请先运行PPI分析")
else:
    print(f"\n找到 {len(ppi_objects)} 个PPI distance对象\n")
    
    # 统计距离分布
    distance_ranges = {
        "0-3Å": 0,
        "3-5Å": 0,
        "5-7Å": 0,
        "7-10Å": 0,
        ">10Å": 0
    }
    
    for obj in sorted(ppi_objects):
        # 获取distance对象的信息
        try:
            # 获取distance值（PyMOL内部存储）
            # 这个方法可能不太准确，但可以给我们一个概念
            info = cmd.get_distance(obj)
            
            if info:
                dist = info
                print(f"{obj}: {dist:.2f} Å")
                
                # 分类
                if dist <= 3.0:
                    distance_ranges["0-3Å"] += 1
                elif dist <= 5.0:
                    distance_ranges["3-5Å"] += 1
                elif dist <= 7.0:
                    distance_ranges["5-7Å"] += 1
                elif dist <= 10.0:
                    distance_ranges["7-10Å"] += 1
                else:
                    distance_ranges[">10Å"] += 1
        except:
            print(f"{obj}: 无法获取距离")
    
    print("\n" + "="*60)
    print("距离分布统计:")
    print("="*60)
    for range_name, count in distance_ranges.items():
        print(f"  {range_name}: {count} 个")
    
    print("\n⚠️ 如果有>5Å的距离，说明可能还在使用CA-CA距离")
    print("   正常的相互作用距离应该都在5Å以内")

print("="*60 + "\n")
