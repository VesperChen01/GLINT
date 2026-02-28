#!/usr/bin/env python3
"""
Comprehensive Chinese to English translation script for GLINT codebase.
Handles technical comments, docstrings, and log messages.
"""

import re
import os
from pathlib import Path

# Comprehensive translation dictionary - organized by category
TRANSLATIONS = {
    # Core technical terms
    "三元复合物": "ternary complex",
    "三元复合体": "ternary complex",
    "分子胶": "molecular glue",
    "蛋白质": "protein",
    "蛋白": "protein",
    "残基": "residue",
    "原子": "atom",
    "配体": "ligand",
    "受体": "receptor",
    "供体": "donor",
    "氢键": "hydrogen bond",
    "盐桥": "salt bridge",
    "疏水": "hydrophobic",
    "亲水": "hydrophilic",
    "相互作用": "interaction",
    "结合": "binding",
    "口袋": "pocket",
    "表面": "surface",
    "表面积": "surface area",
    "埋藏表面积": "buried surface area",
    "坐标": "coordinates",
    "距离": "distance",
    "角度": "angle",
    "能量": "energy",
    "电荷": "charge",
    "极性": "polar",
    "非极性": "nonpolar",
    "芳香": "aromatic",
    "环": "ring",
    "链": "chain",
    "模型": "model",
    "结构": "structure",
    "复合物": "complex",
    "对接": "docking",
    "突变": "mutation",
    "优化": "optimization",
    "筛选": "screening",
    "评分": "scoring",
    "分析": "analysis",
    "检测": "detection",
    "可视化": "visualization",
    "计算": "calculation",
    "评估": "evaluation",
    
    # Analysis terms
    "界面": "interface",
    "接触": "contact",
    "接触数": "contact count",
    "质心": "center of mass",
    "重心": "center of geometry",
    "几何": "geometry",
    "几何中心": "geometric center",
    "特征": "feature",
    "性质": "property",
    "参数": "parameter",
    "指标": "metric",
    "阈值": "threshold",
    "截断": "cutoff",
    
    # Common phrases - conditionals
    "如果": "if",
    "否则": "else",
    "当": "when",
    "则": "then",
    "或": "or",
    "和": "and",
    "但": "but",
    "因为": "because",
    "所以": "therefore",
    "由于": "due to",
    "通过": "through",
    "使用": "using",
    "基于": "based on",
    "根据": "according to",
    "对于": "for",
    "关于": "about",
    "在": "in",
    "从": "from",
    "到": "to",
    "与": "with",
    "为": "as",
    
    # Common verbs
    "计算": "calculate",
    "提取": "extract",
    "获取": "get",
    "返回": "return",
    "生成": "generate",
    "创建": "create",
    "删除": "delete",
    "更新": "update",
    "加载": "load",
    "保存": "save",
    "导入": "import",
    "导出": "export",
    "处理": "process",
    "解析": "parse",
    "转换": "convert",
    "匹配": "match",
    "查找": "find",
    "搜索": "search",
    "过滤": "filter",
    "排序": "sort",
    "合并": "merge",
    "分离": "separate",
    "分组": "group",
    "统计": "count",
    "判断": "determine",
    "检查": "check",
    "验证": "validate",
    "确认": "confirm",
    "显示": "display",
    "隐藏": "hide",
    "高亮": "highlight",
    "标记": "mark",
    "着色": "color",
    "渲染": "render",
    "绘制": "draw",
    "输出": "output",
    "打印": "print",
    "记录": "log",
    "报告": "report",
    "警告": "warning",
    "错误": "error",
    "异常": "exception",
    "失败": "failed",
    "成功": "success",
    "完成": "completed",
    "开始": "start",
    "结束": "end",
    "运行": "run",
    "执行": "execute",
    "调用": "call",
    "应用": "apply",
    "设置": "set",
    "配置": "configure",
    "初始化": "initialize",
    "重置": "reset",
    "清空": "clear",
    "启用": "enable",
    "禁用": "disable",
    "选择": "select",
    "取消": "cancel",
    "跳过": "skip",
    "继续": "continue",
    "停止": "stop",
    "暂停": "pause",
    "恢复": "resume",
}

def translate_text(text):
    """Apply translations to text."""
    # Sort by length (longest first) to avoid partial matches
    for chinese, english in sorted(TRANSLATIONS.items(), key=lambda x: len(x[0]), reverse=True):
        text = text.replace(chinese, english)
    return text

def process_file(filepath):
    """Process a single Python file."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Check if file contains Chinese
        if not re.search(r'[\u4e00-\u9fff]', content):
            return False, "No Chinese found"
        
        original_content = content
        content = translate_text(content)
        
        if content != original_content:
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(content)
            return True, "Translated"
        
        return False, "No changes"
    
    except Exception as e:
        return False, f"Error: {e}"

def main():
    """Main translation process."""
    glint_dir = Path('glint')
    modified_files = []
    skipped_files = []
    
    # Process all Python files in glint directory
    for py_file in glint_dir.rglob('*.py'):
        modified, status = process_file(py_file)
        if modified:
            modified_files.append(str(py_file))
            print(f"✓ {py_file}")
        elif "Chinese" in status:
            skipped_files.append(str(py_file))
    
    print(f"\n{'='*60}")
    print(f"Translation Summary:")
    print(f"  Modified: {len(modified_files)} files")
    print(f"  Skipped (no Chinese): {len(skipped_files)} files")
    print(f"{'='*60}")
    
    # Check for remaining Chinese
    print("\nChecking for remaining Chinese characters...")
    remaining = []
    for py_file in glint_dir.rglob('*.py'):
        try:
            with open(py_file, 'r', encoding='utf-8') as f:
                if re.search(r'[\u4e00-\u9fff]', f.read()):
                    remaining.append(str(py_file))
        except:
            pass
    
    if remaining:
        print(f"\n⚠ {len(remaining)} files still contain Chinese:")
        for f in remaining[:10]:  # Show first 10
            print(f"  - {f}")
        if len(remaining) > 10:
            print(f"  ... and {len(remaining) - 10} more")
    else:
        print("\n✓ All Chinese characters have been translated!")

if __name__ == '__main__':
    main()

