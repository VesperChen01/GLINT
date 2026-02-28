#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Translate Chinese comments and strings to English in Python files
"""

import re
import os
from pathlib import Path

# Translation dictionary for common terms
TRANSLATIONS = {
    # Configuration
    "配置": "Configuration",
    "依赖包": "Dependencies",
    "注意": "Note",
    "已移除": "removed",
    "因为": "because",
    "不可用": "unavailable",
    "与": "with",
    "版本": "version",
    "保持一致": "consistent",
    
    # Common actions
    "查找": "Find",
    "检查": "Check",
    "安装": "Install",
    "路径": "path",
    "常见": "Common",
    "获取": "Get",
    "源码": "source code",
    "目录": "directory",
    "优先": "priority",
    "同级": "sibling",
    "打包后的": "packaged",
    "临时": "temporary",
    
    # GUI related
    "设置": "Set",
    "窗口": "window",
    "图标": "icon",
    "居中": "Center",
    "从": "from",
    "动态": "dynamically",
    "版本": "version",
    "写入": "Write",
    "日志": "log",
    "线程安全": "thread-safe",
    "浏览": "Browse",
    "环境": "environment",
    "开始": "Start",
    "执行": "Execute",
    "创建": "Create",
    "删除": "Delete",
    "旧": "old",
    "复制": "Copy",
    "文件": "files",
    "快捷方式": "shortcut",
    "完成": "Complete",
    
    # Analysis related
    "导出": "Export",
    "打印": "Print",
    "汇总": "Summary",
    "检测": "Detection",
    "总": "Total",
    "结构数": "structures",
    "成功": "Successful",
    "分析": "analysis",
    "数": "count",
    "批量": "Batch",
    "界面": "interface",
    "参数": "Parameters",
    "文件路径": "file path",
    "列表": "list",
    "链对": "chain pairs",
    "如果": "If",
    "只有": "only",
    "一个": "one",
    "元组": "tuple",
    "则": "then",
    "应用于": "apply to",
    "所有": "all",
    "距离": "distance",
    "阈值": "threshold",
    "输出": "output",
    
    # Technical terms
    "静电互补性": "electrostatic complementarity",
    "核心": "core",
    "依赖": "dependency",
    "在": "on",
    "上": "on",
    "不稳定": "unstable",
    "改用": "use instead",
    "移除该": "remove this",
    "验证": "Verify",
    "是否": "whether",
    "成功": "successful",
    "的": "",
}

def translate_line(line):
    """Translate Chinese in a single line"""
    # Skip lines that are mostly code
    if line.strip().startswith(('import ', 'from ', 'class ', 'def ', '@')):
        return line
    
    # Translate common patterns
    for cn, en in TRANSLATIONS.items():
        if cn in line:
            line = line.replace(cn, en)
    
    return line

def process_file(filepath):
    """Process a single Python file"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Check if file has Chinese
        has_chinese = any(re.search(r'[\u4e00-\u9fff]', line) for line in lines)
        if not has_chinese:
            return False
        
        # Translate lines
        translated_lines = [translate_line(line) for line in lines]
        
        # Write back
        with open(filepath, 'w', encoding='utf-8') as f:
            f.writelines(translated_lines)
        
        return True
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return False

def main():
    """Main function"""
    glint_dir = Path('glint')
    
    if not glint_dir.exists():
        print("Error: glint directory not found")
        return
    
    # Find all Python files
    py_files = list(glint_dir.rglob('*.py'))
    
    print(f"Found {len(py_files)} Python files")
    print("Processing...")
    
    processed = 0
    for py_file in py_files:
        if process_file(py_file):
            processed += 1
            print(f"  Processed: {py_file}")
    
    print(f"\nCompleted: {processed} files translated")

if __name__ == '__main__':
    main()

