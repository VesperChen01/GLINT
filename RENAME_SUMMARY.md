# 🎉 GlueTK Renaming Summary

## ✅ 重命名完成

插件已成功从 **MolStruct** 重命名为 **GlueTK** (Glue ToolKit)。

---

## 📝 主要变更

### 1. 品牌名称
- **旧名称**: MolStruct Plugin for PyMOL
- **新名称**: GlueTK - PyMOL Plugin for Molecular Glue Analysis
- **定位**: Molecular Glue vs PROTAC Classification Toolkit

### 2. GUI 启动命令
- **新命令**: `gluetk_gui` ⭐ (推荐使用)
- **旧命令**: `molstruct_gui` (保留为兼容性别名)

### 3. PyMOL 菜单
- **新菜单**: Plugins → **GlueTK - Molecular Glue Analyzer**

### 4. 主要类名
- `MolStructDialog` → `GlueTKDialog`

---

## 📂 修改的文件

### 核心文件
- ✅ `gluetk/__init__.py` - 插件入口，命令注册
- ✅ `gluetk/unified_gui.py` - GUI 主窗口
- ✅ `gluetk/README.md` - 插件文档

### 文档文件
- ✅ `WARP.md` - 开发文档
- ✅ `manuscript/JCIM_manuscript_outline.md` - 论文大纲
- ✅ `manuscript/figure_generation_guide.md` - 图表生成指南
- ✅ `manuscript/submission_checklist.md` - 投稿清单

### 脚本文件
- ✅ `validation/benchmark_analysis.py` - 验证脚本
- ✅ `reload_plugin.py` - 重载脚本
- ✅ `gluetk/*.py` - 所有 Python 模块
- ✅ `gluetk/*.sh` - Shell 脚本

### 测试文件
- ✅ `test_heatmap_example/README.md`

---

## 🔄 兼容性保留

为了向后兼容，以下内容**保持不变**：

### 1. 目录名称
- ✅ `gluetk/` - **已重命名** (原 `molstruct_plugin/`)
- 所有导入路径已更新为 `from gluetk import ...`

### 2. 命令别名
- ✅ `molstruct_gui` - 作为 `gluetk_gui` 的别名继续可用
- 用户可以使用任一命令启动 GUI

---

## 🚀 使用方法

### 在 PyMOL 中加载插件

```python
# 方法 1: 通过 PyMOL Plugin Manager
# Plugin → Plugin Manager → Install New Plugin → 选择 gluetk 目录

# 方法 2: 手动加载
run /path/to/glue-pymol/gluetk/__init__.py

# 例如：
run /Users/vesper/Desktop/git/glue-pymol/gluetk/__init__.py
```

### 启动 GUI

```python
# 推荐：使用新命令
gluetk_gui

# 或使用兼容命令
molstruct_gui
```

### 菜单访问
```
Plugins → GlueTK - Molecular Glue Analyzer
```

---

## 📊 功能完整性

所有原有功能**完全保留**，包括：

### ✨ 分子胶专用分析
- 🔗 PPI Interface Analysis
- ✨ Neo-Epitope Detection
- 🎯 Molecular Glue vs PROTAC Classification
- 🧬 G-Motif Recognition
- 🔬 G-Motif Glue Binding Validation

### 🔬 通用相互作用分析
- Hydrogen Bonds (≤2.8Å Schrödinger standard)
- Salt Bridges (≤4.0Å)
- Hydrophobic Interactions
- π-π Stacking, π-Cation
- Metal Coordination

### ⚡ 结合能评分
- Binary Complex Scoring
- Ternary Complex Scoring (with cooperativity)
- Binding Heatmap Generation

### 🎨 可视化
- 3D PyMOL Visualization
- 2D Interaction Diagrams
- Network Plots
- Electrostatics (APBS/Quick)

---

## ⚠️ 重要提示

### 对于现有用户
1. **需要更新导入路径** - `from molstruct_plugin import ...` → `from gluetk import ...`
2. **命令保持兼容** - `molstruct_gui` 仍然可用
3. **建议全面迁移** - 使用新名称 `gluetk_gui` 和 `from gluetk import ...`

### 对于新用户
1. **推荐使用** `gluetk_gui` 命令
2. **查看文档** 时认准 **GlueTK** 品牌
3. **引用论文** 时使用 GlueTK 名称

---

## 📚 GitHub 仓库

- **新地址**: https://github.com/VesperChen01/GlueTK.git
- **Star ⭐ us**: 如果 GlueTK 对你的研究有帮助！

---

## 🔮 未来计划

### 已完成优化
1. **目录重命名**: ✅ 完成
   - `molstruct_plugin/` → `gluetk/`
   - 所有导入路径已更新
   - 环境名称已更新为 `gluetk`

2. **命令简化**: 可考虑添加更短的别名
   - `gtk_gui` 作为 `gluetk_gui` 的超短别名
   - `gtk.analyze_ppi()` 等命名空间组织

---

## ✅ 验证清单

在 PyMOL 中测试以下功能：

- [ ] `gluetk_gui` 命令成功启动 GUI
- [ ] `molstruct_gui` 兼容命令正常工作
- [ ] GUI 窗口标题显示 "GlueTK - Molecular Glue Analyzer"
- [ ] About 页面显示 GlueTK 品牌
- [ ] 所有分析功能正常运行
- [x] 导入语句 `from gluetk import ...` 正常工作
- [ ] PyMOL 菜单显示 "GlueTK - Molecular Glue Analyzer"

---

## 📞 联系方式

- **Developer**: Roufen Chen
- **Email**: 12319021@zju.edu.cn
- **GitHub**: https://github.com/VesperChen01/GlueTK

---

**版本**: 1.0.0  
**重命名日期**: 2025-11-09  
**状态**: ✅ 完成
