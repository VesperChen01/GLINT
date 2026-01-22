# GLINT 修复日志

## 2024年修复记录

### 1. 移除 PPI 可视化中的中文图例 ✅

**问题描述：**
在蛋白质-蛋白质相互作用（PPI）分析后，PyMOL 可视化界面中会显示中文的相互作用类型图例伪原子，包括：
- 氢键: X
- 盐桥: X
- 疏水接触: X
- π-π堆积: X
- 阳离子-π: X

**修复内容：**
- 文件：`glint/ppi_analyzer.py`
- 修改位置：第 908-919 行
- 删除了创建图例伪原子（pseudoatom）的代码
- 保留了终端统计信息输出
- 删除了不再需要的 `TYPE_LABEL_MAP` 映射

**修改后效果：**
- ✅ PyMOL 可视化界面中不再显示任何图例文字
- ✅ 相互作用的虚线和颜色仍然正常显示
- ✅ 统计信息仍然在终端输出，方便查看
- ✅ 界面更加简洁专业

---

### 2. 修复 PPI 分析不显示残基标签的问题 ✅

**问题描述：**
PPI 相互作用分析后，高亮的残基不显示标签（如 A123, B456 等）。

**原因分析：**
代码中使用 `cmd.label()` 设置了标签文本，但没有显式调用 `cmd.show("labels")` 来显示标签。在某些 PyMOL 配置下，标签默认是隐藏的。

**修复内容：**
- 文件：`glint/ppi_analyzer.py`
- 修改位置：第 691 行
- 添加了 `cmd.show("labels", ca_sel)` 来显式显示残基标签

**修改前：**
```python
cmd.label(ca_sel, f"'{label_text}'")
cmd.set("label_color", "black", ca_sel)
cmd.set("label_size", 14, ca_sel)
```

**修改后：**
```python
cmd.label(ca_sel, f"'{label_text}'")
cmd.set("label_color", "black", ca_sel)
cmd.set("label_size", 14, ca_sel)
cmd.show("labels", ca_sel)  # 显式显示标签
```

**修改后效果：**
- ✅ PPI 分析后，相互作用残基会显示标签（如 A123, B456）
- ✅ 标签显示为黑色，字体大小为 14
- ✅ 标签位置在 CA 原子上
- ✅ 使用单字母氨基酸代码 + 残基号（如 A123 表示 ALA 123）

---

### 3. 改进 PyMOL 安装策略，避免依赖冲突 ✅

**问题描述：**
当用户系统中已有 PyMOL.app 时，安装脚本会直接使用它，但这会导致依赖冲突：
- 系统 PyMOL 使用自己的 Python 环境
- GLINT 的依赖安装在 conda 环境中
- 两者无法互通，导致 GLINT 无法正常工作

**修复内容：**
- 文件：`install_glint.sh`
- 修改了 PyMOL 检测和安装逻辑（第 66-122 行）
- 修改了启动器生成逻辑（第 351-407 行）
- 更新了安装完成提示信息（第 413-439 行）

**新的安装策略：**

1. **优先使用 conda 环境中的 PyMOL**（推荐）
   - 即使系统已有 PyMOL.app，也会在 conda 环境中安装独立的 PyMOL
   - 所有依赖都在同一个 conda 环境中，完全隔离
   - 避免与系统 PyMOL 冲突
   - 用户的系统 PyMOL 不受影响

2. **备用方案：使用系统 PyMOL.app**
   - 仅在 conda 安装 PyMOL 失败时使用
   - 会尝试将 conda 环境的 site-packages 添加到 PYTHONPATH
   - 显示警告信息，提示可能存在依赖冲突

**启动器改进：**

- **Conda PyMOL 启动器**（推荐）：
  ```bash
  conda activate glint
  pymol -d "import glint; glint.glint_gui()"
  ```

- **系统 PyMOL.app 启动器**（备用）：
  ```bash
  export PYTHONPATH="$CONDA_ENV_PYTHON:$PYTHONPATH"
  /Applications/PyMOL.app/Contents/MacOS/PyMOL -d "
  import sys
  sys.path.insert(0, conda_site_packages)
  import glint
  glint.glint_gui()
  "
  ```

**用户体验改进：**

安装时的提示信息：
```
✅ 检测到系统 PyMOL.app: /Applications/PyMOL.app/Contents/MacOS/PyMOL
   在 conda 环境中安装 PyMOL（推荐方式）...
✅ PyMOL 已安装到 conda 环境
💡 提示: 检测到系统已有 PyMOL.app，但 GLINT 将使用 conda 环境中的 PyMOL
   这样可以避免依赖冲突，确保所有功能正常工作
```

安装完成后的提示：
```
🎉 GLINT 安装成功！

PyMOL 类型: conda
✅ 使用 conda 环境中的 PyMOL（推荐）
   所有依赖已隔离，避免与系统 PyMOL 冲突
💡 提示: 检测到系统 PyMOL.app，但 GLINT 使用独立的 conda PyMOL
   这样可以确保依赖兼容性，您的系统 PyMOL 不会受影响
```

**优势：**
- ✅ 完全避免依赖冲突
- ✅ 用户的系统 PyMOL 不受影响
- ✅ 所有 GLINT 功能都能正常工作
- ✅ 环境隔离，便于维护和调试
- ✅ 即使系统 PyMOL 更新也不会影响 GLINT

---

## 测试建议

### 测试场景 1：全新安装（无系统 PyMOL）
```bash
bash install_glint.sh
```
预期结果：
- 在 conda 环境中安装 PyMOL
- 创建 GLINT.app 启动器
- 双击启动器可正常打开 GLINT

### 测试场景 2：已有系统 PyMOL.app
```bash
# 确保 /Applications/PyMOL.app 存在
bash install_glint.sh
```
预期结果：
- 检测到系统 PyMOL.app
- 仍然在 conda 环境中安装独立的 PyMOL
- 显示提示信息说明使用 conda PyMOL
- 双击启动器可正常打开 GLINT
- 系统 PyMOL.app 不受影响

### 测试场景 3：conda 安装 PyMOL 失败
```bash
# 模拟 conda 安装失败
bash install_glint.sh
```
预期结果：
- 如果有系统 PyMOL.app，回退使用它
- 显示警告信息
- 尝试注入 conda 环境路径
- 可能存在部分功能不可用

---

## 文件修改清单

1. `glint/ppi_analyzer.py` - 移除中文图例 + 修复残基标签显示
2. `install_glint.sh` - 改进 PyMOL 安装策略
3. `sync_to_pymol.sh` - 新增：同步修改到 PyMOL 启动目录的便捷脚本

## 如何应用修改

### 方法 1：使用同步脚本（推荐）
```bash
# 在项目根目录运行
bash sync_to_pymol.sh
```

这会自动将所有修改同步到 `~/.pymol/startup/glint/` 目录。

### 方法 2：重新安装
```bash
bash install_glint.sh
```

### 方法 3：手动复制
```bash
cp glint/ppi_analyzer.py ~/.pymol/startup/glint/ppi_analyzer.py
```

### 应用修改后
1. **如果 PyMOL 正在运行**：
   - 重启 PyMOL
   - 或在 PyMOL 命令行运行：`reinitialize`
   - 然后重新导入：`import glint; glint.glint_gui()`

2. **如果 PyMOL 未运行**：
   - 直接启动 PyMOL 即可使用更新后的版本

## 验证修改是否生效

运行 PPI 分析后，检查：
- ✅ 残基标签正确显示（如 A123, B456）
- ✅ 不显示中文图例（氢键、盐桥等）
- ✅ 相互作用虚线正常显示

## 兼容性

- ✅ macOS (主要测试平台)
- ⚠️ Linux (需要测试)
- ⚠️ Windows (需要测试)

## 后续建议

1. 考虑为 Linux 和 Windows 创建类似的安装脚本改进
2. 添加环境检测工具，帮助用户诊断依赖问题
3. 考虑提供"修复模式"，重新安装 conda PyMOL

