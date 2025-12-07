"""
列出所有PyMOL对象

在PyMOL中运行:
run /Users/vesper/Desktop/git/GlueTK/GlueTK/list_all_objects.py
"""

from pymol import cmd

print("\n" + "="*60)
print("所有PyMOL对象")
print("="*60)

all_objects = cmd.get_names("all")
print(f"\n总共 {len(all_objects)} 个对象:\n")

# 按前缀分组
groups = {}
for obj in all_objects:
    prefix = obj.split("_")[0] if "_" in obj else obj
    if prefix not in groups:
        groups[prefix] = []
    groups[prefix].append(obj)

for prefix in sorted(groups.keys()):
    objs = groups[prefix]
    print(f"\n{prefix}* : {len(objs)} 个对象")
    for obj in objs[:5]:  # 只显示前5个
        print(f"  - {obj}")
    if len(objs) > 5:
        print(f"  ... 还有 {len(objs)-5} 个")

print("\n" + "="*60 + "\n")
