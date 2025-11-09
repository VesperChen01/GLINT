# 🔄 GlueTK Migration Guide

## 从 MolStruct 迁移到 GlueTK

如果你之前使用过 MolStruct，本指南将帮助你快速迁移到 GlueTK。

---

## 📝 主要变更

### 1. 目录名称变更
```bash
# 旧目录
molstruct_plugin/

# 新目录
gluetk/
```

### 2. 导入路径变更
```python
# ❌ 旧的导入方式（不再有效）
from molstruct_plugin.interaction_analyzer import analyze_pdb_interactions
from molstruct_plugin import __version__

# ✅ 新的导入方式
from gluetk.interaction_analyzer import analyze_pdb_interactions
from gluetk import __version__
```

### 3. 命令变更
```python
# ✅ 推荐：使用新命令
gluetk_gui

# ✅ 兼容：旧命令仍可用
molstruct_gui  # 会自动调用 gluetk_gui
```

---

## 🚀 快速迁移步骤

### Step 1: 更新 PyMOL 插件路径

如果你之前手动加载插件：

```python
# ❌ 旧的加载方式
run /path/to/molstruct_plugin/__init__.py

# ✅ 新的加载方式
run /path/to/gluetk/__init__.py
```

### Step 2: 更新 Python 脚本中的导入

使用查找替换功能批量更新：

```bash
# 在你的项目目录中运行
find . -name "*.py" -type f -exec sed -i '' 's/from molstruct_plugin/from gluetk/g' {} \;
```

或者手动更新每个文件：

```python
# 修改前
from molstruct_plugin.ppi_analyzer import analyze_protein_protein_interface
from molstruct_plugin.binding_score import calculate_ternary_score

# 修改后
from gluetk.ppi_analyzer import analyze_protein_protein_interface
from gluetk.binding_score import calculate_ternary_score
```

### Step 3: 更新配置文件和脚本

如果你的脚本中包含硬编码的路径：

```python
# 修改前
plugin_dir = "/path/to/molstruct_plugin"
sys.path.insert(0, "/path/to/molstruct_plugin")

# 修改后
plugin_dir = "/path/to/gluetk"
sys.path.insert(0, "/path/to/gluetk")
```

---

## ✅ 验证迁移

### 测试 1: 验证导入
```python
import sys
sys.path.insert(0, '/path/to/gluetk')

from gluetk import __version__, __author__
print(f"GlueTK v{__version__} by {__author__}")
# 应该输出: GlueTK v1.0.0 by Vesper
```

### 测试 2: 验证 PyMOL 命令
在 PyMOL 中运行：
```python
# 加载插件
run /path/to/gluetk/__init__.py

# 测试新命令
gluetk_gui

# 测试兼容命令
molstruct_gui
```

### 测试 3: 验证功能模块
```python
from gluetk.interaction_analyzer import analyze_pdb_interactions
from gluetk.ppi_analyzer import analyze_protein_protein_interface
from gluetk.binding_score import calculate_ternary_score

# 如果没有报错，说明迁移成功！
```

---

## 🔍 常见问题

### Q1: 我能同时使用 molstruct_plugin 和 gluetk 吗？
**A**: 不能。目录已经重命名，`molstruct_plugin` 不再存在。请完全迁移到 `gluetk`。

### Q2: molstruct_gui 命令还能用吗？
**A**: 可以！`molstruct_gui` 作为兼容别名保留，它会自动调用 `gluetk_gui`。但建议使用新命令。

### Q3: 我的旧脚本会崩溃吗？
**A**: 如果你的脚本中使用了 `from molstruct_plugin import ...`，需要更新为 `from gluetk import ...`。

### Q4: Conda 环境名称改了吗？
**A**: 是的，环境名称从 `molstruct` 改为 `gluetk`。如果你使用自定义环境，需要更新：
```bash
# 旧环境
conda activate molstruct

# 新环境（如果使用 env_checker 自动创建）
conda activate gluetk
```

### Q5: 我需要重新安装依赖吗？
**A**: 不需要。如果你已经安装了 rdkit, scipy, matplotlib, numpy 等依赖，它们仍然可用。

---

## 📚 更新的文档路径

### 插件文档
- `gluetk/README.md` - 主要文档
- `gluetk/QUICK_START.md` - 快速开始

### 开发文档
- `WARP.md` - 开发指南
- `RENAME_SUMMARY.md` - 重命名总结
- `MIGRATION_GUIDE.md` - 本迁移指南

### 论文相关
- `manuscript/JCIM_manuscript_outline.md`
- `manuscript/figure_generation_guide.md`
- `manuscript/submission_checklist.md`

---

## 🆘 迁移遇到问题？

### 错误: ModuleNotFoundError: No module named 'molstruct_plugin'

**解决方案**:
```python
# 检查导入路径
# ❌ 错误
from molstruct_plugin import something

# ✅ 正确
from gluetk import something
```

### 错误: FileNotFoundError: /path/to/molstruct_plugin/__init__.py

**解决方案**:
```python
# 更新加载路径
# ❌ 错误
run /path/to/molstruct_plugin/__init__.py

# ✅ 正确
run /path/to/gluetk/__init__.py
```

### GUI 无法启动

**解决方案**:
```python
# 1. 清除旧模块缓存
import sys
modules_to_remove = [k for k in sys.modules.keys() 
                     if 'molstruct' in k.lower() or 'gluetk' in k.lower()]
for mod in modules_to_remove:
    del sys.modules[mod]

# 2. 重新加载
run /path/to/gluetk/__init__.py

# 3. 启动 GUI
gluetk_gui
```

---

## 📞 需要帮助？

- **GitHub Issues**: https://github.com/VesperChen01/GlueTK/issues
- **Email**: 12319021@zju.edu.cn

---

**最后更新**: 2025-11-09  
**适用版本**: GlueTK v1.0.0
