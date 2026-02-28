#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to replace all Chinese text with English in Python files
"""

import re
import os
import sys

# Comprehensive translation dictionary
TRANSLATIONS = {
    # Module descriptions
    "三元复合物评估模块": "Ternary Complex Evaluator Module",
    "整合三类计算性质": "Integrates three types of computational properties",
    "表面积计算": "surface area calculation",
    "接触数统计": "contact count statistics",
    "小分子理化参数": "small molecule physicochemical parameters",
    "三元复合物重心": "ternary complex center of mass",
    "几何特征": "geometric features",
    "依赖": "Dependencies",
    
    # Data structures
    "数据结构": "Data Structures",
    "原子信息": "Atom information",
    "链信息": "Chain information",
    "计算质心": "Calculate center of mass",
    "界面特征": "Interface features",
    "总掩埋表面积": "Total buried surface area",
    "界面": "interface",
    "接触数": "contact count",
    "配体特征": "Ligand features",
    "距离特征": "Distance features",
    "最小原子距离": "Minimum atomic distance",
    "偏移": "shift",
    "连线的垂直距离": "perpendicular distance to line",
    "夹角": "angle",
    "度": "degrees",
    "三元复合物综合特征": "Comprehensive ternary complex features",
    "双面性指数": "Duality index",
    "计算器": "Calculator",
    "使用": "using",
    
    # Common actions
    "导出": "Export",
    "导入": "Import",
    "导航栏页面索引常量": "Navigation bar page index constants",
    "消除魔法数字": "Eliminate magic numbers",
    "加载": "Load",
    "保存": "Save",
    "删除": "Delete",
    "更新": "Update",
    "创建": "Create",
    "初始化": "Initialize",
    "配置": "Configuration",
    "设置": "Settings",
    "选项": "Options",
    "参数": "Parameters",
    "结果": "Results",
    "错误": "Error",
    "警告": "Warning",
    "信息": "Information",
    "成功": "Success",
    "失败": "Failed",
    "完成": "Completed",
    "取消": "Cancel",
    "确认": "Confirm",
    "关闭": "Close",
    "打开": "Open",
    "开始": "Start",
    "停止": "Stop",
    "暂停": "Pause",
    "继续": "Continue",
    "重试": "Retry",
    "跳过": "Skip",
    "返回": "Return",
    "退出": "Exit",
    "帮助": "Help",
    "关于": "About",
    "版本": "Version",
    "作者": "Author",
    "日期": "Date",
    "时间": "Time",
    "文件": "File",
    "目录": "Directory",
    "路径": "Path",
    "名称": "Name",
    "类型": "Type",
    "大小": "Size",
    "状态": "Status",
    "进度": "Progress",
    "描述": "Description",
    "备注": "Note",
    "标签": "Label",
    "分类": "Category",
    "搜索": "Search",
    "过滤": "Filter",
    "排序": "Sort",
    "显示": "Display",
    "隐藏": "Hide",
    "启用": "Enable",
    "禁用": "Disable",
    "选择": "Select",
    "全选": "Select All",
    "清空": "Clear",
    "重置": "Reset",
    "应用": "Apply",
    "预览": "Preview",
    "编辑": "Edit",
    "查看": "View",
    "复制": "Copy",
    "粘贴": "Paste",
    "剪切": "Cut",
    "撤销": "Undo",
    "重做": "Redo",
    "刷新": "Refresh",
    "下载": "Download",
    "上传": "Upload",
    "分享": "Share",
    "打印": "Print",
    "导出为": "Export as",
    "导入自": "Import from",
    "另存为": "Save as",
    "新建": "New",
    "添加": "Add",
    "移除": "Remove",
    "修改": "Modify",
    "替换": "Replace",
    "查找": "Find",
    "定位": "Locate",
    "跳转": "Jump to",
    "上一个": "Previous",
    "下一个": "Next",
    "第一个": "First",
    "最后一个": "Last",
    "总计": "Total",
    "数量": "Count",
    "序号": "Index",
    "编号": "Number",
    "标识": "ID",
    "键": "Key",
    "值": "Value",
    "属性": "Property",
    "方法": "Method",
    "函数": "Function",
    "类": "Class",
    "模块": "Module",
    "包": "Package",
    "库": "Library",
    "依赖项": "Dependency",
    "组件": "Component",
    "插件": "Plugin",
    "扩展": "Extension",
    "工具": "Tool",
    "实用程序": "Utility",
}

def replace_chinese_in_file(filepath):
    """Replace Chinese text with English in a single file"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Check if file has Chinese
        if not re.search(r'[\u4e00-\u9fff]', content):
            return False
        
        original_content = content
        
        # Apply translations
        for chinese, english in TRANSLATIONS.items():
            content = content.replace(chinese, english)
        
        # Write back if changed
        if content != original_content:
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(content)
            return True
        
        return False
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return False

def main():
    """Main function to process all Python files"""
    processed_count = 0
    modified_count = 0
    
    # Find all Python files
    for root, dirs, files in os.walk('glint'):
        for file in files:
            if file.endswith('.py'):
                filepath = os.path.join(root, file)
                processed_count += 1
                
                if replace_chinese_in_file(filepath):
                    modified_count += 1
                    print(f"✓ Modified: {filepath}")
    
    print(f"\n{'='*60}")
    print(f"Processed {processed_count} Python files")
    print(f"Modified {modified_count} files")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()

