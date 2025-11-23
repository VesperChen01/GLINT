#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
调试原子连接线显示问题
Debug atom interaction lines
"""

from pymol import cmd
import csv
import os

def debug_atom_lines():
    """调试为什么原子连接线不显示"""
    
    print("\n" + "="*60)
    print("调试原子连接线 (Debug Atom Lines)")
    print("="*60 + "\n")
    
    # 1. 检查CSV文件内容
    csv_file = "test_interactions.csv"
    if not os.path.exists(csv_file):
        print(f"❌ CSV文件不存在: {csv_file}")
        print("   请先运行: analyze_pdb_interactions('your_structure', output_csv='test_interactions.csv')")
        return
    
    print("步骤 1: 检查CSV文件内容...")
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        
    print(f"CSV 包含 {len(rows)} 行数据")
    
    # 检查列名
    if rows:
        headers = list(rows[0].keys())
        print(f"CSV 列名: {headers}")
        
        # 检查是否有原子信息
        has_atom1 = any('atom' in k.lower() for k in headers)
        has_atom2 = any('atom2' in k.lower() for k in headers)
        
        if has_atom1 and has_atom2:
            print("✅ CSV包含原子信息列")
        else:
            print("❌ CSV缺少原子信息列")
    
    # 2. 显示前3行的原子信息
    print("\n步骤 2: 检查原子信息内容...")
    for i, row in enumerate(rows[:3], 1):
        atom1 = None
        atom2 = None
        
        # 尝试多种列名
        for key in ['Atom1', 'atom1', 'Ligand_Atom', 'ligand_atom']:
            if key in row and row[key]:
                atom1 = row[key]
                break
        
        for key in ['Atom2', 'atom2', 'Protein_Atom', 'protein_atom']:
            if key in row and row[key]:
                atom2 = row[key]
                break
        
        print(f"Row {i}:")
        print(f"  Atom1: '{atom1}' (empty: {not atom1})")
        print(f"  Atom2: '{atom2}' (empty: {not atom2})")
        
        # 显示完整行数据
        if i == 1:
            print("  完整数据:", dict(row))
    
    # 3. 手动测试创建距离对象
    print("\n步骤 3: 测试创建距离对象...")
    
    # 获取当前对象
    objs = cmd.get_names("objects")
    if not objs:
        print("❌ 没有加载的对象")
        return
    
    obj = objs[0]
    print(f"使用对象: {obj}")
    
    # 尝试创建一个测试距离对象
    if rows and len(rows) > 0:
        row = rows[0]
        
        # 解析第一个相互作用
        chain1 = row.get('Chain1', '')
        res1 = row.get('Residue1', '').split()
        chain2 = row.get('Chain2', '')
        res2 = row.get('Residue2', '').split()
        
        if len(res1) >= 2 and len(res2) >= 2:
            resid1 = res1[1]
            resid2 = res2[1]
            
            # 构建选择
            sel1 = f"{obj} and chain {chain1} and resi {resid1}"
            sel2 = f"{obj} and chain {chain2} and resi {resid2}"
            
            # 检查原子数
            n1 = cmd.count_atoms(sel1)
            n2 = cmd.count_atoms(sel2)
            
            print(f"\n测试选择:")
            print(f"  sel1: {sel1} -> {n1} 原子")
            print(f"  sel2: {sel2} -> {n2} 原子")
            
            if n1 > 0 and n2 > 0:
                # 尝试创建距离对象
                try:
                    cmd.distance("test_dist", sel1, sel2)
                    print("✅ 成功创建测试距离对象 'test_dist'")
                    
                    # 设置样式
                    cmd.color("red", "test_dist")
                    cmd.set("dash_width", 5, "test_dist")
                    cmd.show("dashes", "test_dist")
                    
                except Exception as e:
                    print(f"❌ 创建距离对象失败: {e}")
            else:
                print("❌ 选择没有找到原子")
    
    # 4. 检查现有的距离对象
    print("\n步骤 4: 检查现有距离对象...")
    dist_objs = [name for name in cmd.get_names("objects") if name.startswith("dist_")]
    
    if dist_objs:
        print(f"找到 {len(dist_objs)} 个距离对象:")
        for obj_name in dist_objs[:5]:
            print(f"  - {obj_name}")
            # 检查是否可见
            try:
                # PyMOL没有直接的API检查可见性，但可以尝试设置属性
                cmd.show("dashes", obj_name)
                cmd.set("dash_width", 3, obj_name)
            except:
                pass
    else:
        print("❌ 没有找到距离对象")
    
    # 5. 强制显示所有虚线
    print("\n步骤 5: 强制显示所有虚线...")
    cmd.show("dashes")
    cmd.set("dash_width", 3)
    cmd.set("dash_gap", 0.3)
    cmd.set("dash_length", 0.5)
    
    # 列出所有测量对象
    all_dists = cmd.get_names("objects", type=4)  # type=4 是距离对象
    if all_dists:
        print(f"所有距离/测量对象: {all_dists}")
    
    print("\n" + "="*60)
    print("调试完成！")
    print("\n可能的解决方案:")
    print("1. 确保CSV文件包含Atom1和Atom2列")
    print("2. 运行: show dashes")
    print("3. 运行: set dash_width, 3")
    print("4. 检查是否有 'test_dist' 对象（红色粗线）")
    print("="*60 + "\n")

if __name__ == "__main__":
    debug_atom_lines()