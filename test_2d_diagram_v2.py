#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
快速测试 2D 互作图 v2.2 优化效果

使用方法:
1. 在PyMOL中加载蛋白-配体复合物
2. 生成互作CSV文件
3. 运行此脚本:
   python test_2d_diagram_v2.py --csv interactions.csv --ligand MOL

优化亮点:
✅ 圆形氨基酸标签 + 单字母命名
✅ 氨基酸分类着色 (非极性/极性/带电)
✅ 分子完全无氢显示
✅ 疏水作用仅显示半透明圆圈 (无连接线)
✅ 紧凑清爽布局
"""

import argparse
import sys
import os

def test_2d_diagram():
    parser = argparse.ArgumentParser(
        description='测试优化后的2D互作图生成',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--csv', required=True, help='互作CSV文件路径')
    parser.add_argument('--ligand', required=True, help='配体残基名称')
    parser.add_argument('--obj', default=None, help='PyMOL对象名称 (可选)')
    parser.add_argument('--output', default=None, help='输出PNG路径 (可选)')
    parser.add_argument('--protein', default=None, help='蛋白质名称 (用于标题)')
    parser.add_argument('--width', type=int, default=2400, help='图片宽度')
    parser.add_argument('--height', type=int, default=2000, help='图片高度')
    parser.add_argument('--dpi', type=int, default=300, help='分辨率')
    
    args = parser.parse_args()
    
    # 检查CSV文件
    if not os.path.exists(args.csv):
        print(f"❌ 错误: CSV文件不存在: {args.csv}")
        sys.exit(1)
    
    # 导入GlueTK
    try:
        from gluetk.interaction_2d_plot import generate_2d_interaction_diagram
        print("✅ 成功加载 GlueTK 2D互作图模块")
    except ImportError as e:
        print(f"❌ 无法导入GlueTK: {e}")
        print("请确保在GlueTK项目目录运行,或已安装gluetk包")
        sys.exit(1)
    
    # 生成图表
    print(f"\n🎨 开始生成2D互作图...")
    print(f"   配体: {args.ligand}")
    print(f"   CSV: {args.csv}")
    print(f"   分辨率: {args.width}x{args.height} @ {args.dpi} DPI")
    
    output_path = generate_2d_interaction_diagram(
        csv_path=args.csv,
        ligand_resname=args.ligand,
        obj_name=args.obj,
        output_path=args.output,
        width=args.width,
        height=args.height,
        dpi=args.dpi,
        protein_name=args.protein,
        compact=True
    )
    
    if output_path and os.path.exists(output_path):
        print(f"\n✅ 成功! 图片已保存至: {output_path}")
        print("\n🎯 优化特性:")
        print("   • 圆形氨基酸标签 (单字母命名)")
        print("   • 氨基酸分类自动着色")
        print("   • 分子结构完全无氢")
        print("   • 疏水作用无连接线 (仅高亮圆)")
        print("   • 紧凑清爽专业布局")
    else:
        print("\n❌ 生成失败!")
        sys.exit(1)

if __name__ == "__main__":
    test_2d_diagram()
