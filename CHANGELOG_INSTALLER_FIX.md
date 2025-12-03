# 安装器修复日志

## 修复日期
2025-12-03

## 问题描述
GlueTK Installer.app 生成的桌面启动器无法启动，报错：
```
pymol: command not found
```

## 根本原因
1. **conda PyMOL 安装不可靠**：`pymol-open-source` 在 macOS 上通过 conda 安装经常失败
2. **缺少验证**：安装器没有验证 PyMOL 是否真的安装成功
3. **启动器逻辑单一**：只生成一种启动器，假设 PyMOL 必定在 conda 环境中
4. **忽略系统 PyMOL**：没有检测和利用系统已安装的 PyMOL.app

## 修复内容

### 1. 新增独立安装脚本

#### `install_gluetk.sh` - 完整安装脚本
**功能：**
- ✅ 智能检测 PyMOL 安装方式（系统 > conda）
- ✅ 验证 PyMOL 安装是否成功
- ✅ 根据 PyMOL 类型生成对应启动器
- ✅ 清晰的错误提示和建议

**使用场景：** 首次安装或完整重新安装

#### `fix_pymol.sh` - 快速修复脚本
**功能：**
- ✅ 不重装依赖，只修复 PyMOL 问题
- ✅ 优先使用系统 PyMOL.app
- ✅ 必要时在 conda 中安装 PyMOL
- ✅ 重建正确的桌面启动器

**使用场景：** 已运行安装器但无法启动

### 2. 修复 `rebuild_installer.sh`

#### 改进点 1：PyMOL 安装验证
```python
# 第 376-380 行
# 验证 PyMOL 是否安装成功
pymol_check = subprocess.run([self.conda_exe, "run", "-n", ENV_NAME, "which", "pymol"], 
                            capture_output=True, text=True)
if pymol_check.returncode != 0:
    self._log("  ⚠️ PyMOL not installed in conda, will check system PyMOL.app")
```

#### 改进点 2：智能启动器生成
```python
# 第 476-528 行
# 检测 PyMOL 安装方式
pymol_app_path = "/Applications/PyMOL.app/Contents/MacOS/PyMOL"
has_pymol_app = os.path.exists(pymol_app_path) and os.access(pymol_app_path, os.X_OK)

# 检查 conda 环境中是否有 pymol
has_conda_pymol = False
try:
    check_result = subprocess.run([self.conda_exe, "run", "-n", ENV_NAME, "which", "pymol"],
                                capture_output=True, text=True, timeout=5)
    has_conda_pymol = (check_result.returncode == 0)
except:
    pass

# 根据检测结果生成不同的启动器
if has_pymol_app:
    # 使用系统 PyMOL.app + conda 依赖
elif has_conda_pymol:
    # 使用 conda PyMOL
else:
    # 创建错误提示启动器
```

#### 改进点 3：三种启动器模式

**模式 1：系统 PyMOL.app（推荐）**
```bash
#!/bin/bash
# 激活 conda 环境以加载依赖
export CONDA_EXE="/path/to/conda"
eval "$($CONDA_EXE shell.bash hook)" 2>/dev/null
conda activate gluetk 2>/dev/null || true

# 使用系统 PyMOL
"/Applications/PyMOL.app/Contents/MacOS/PyMOL" -d "import sys, os; sys.path.insert(0, os.path.expanduser('~/.pymol/startup')); import gluetk; gluetk.gluetk_gui()"
```

**优点：**
- 性能更好
- 更稳定
- 同时利用 conda 依赖包

**模式 2：conda PyMOL（备选）**
```bash
#!/bin/bash
CONDA_EXE="/path/to/conda"
"$CONDA_EXE" run -n gluetk pymol -d "import sys, os; sys.path.insert(0, os.path.expanduser('~/.pymol/startup')); import gluetk; gluetk.gluetk_gui()"
```

**优点：**
- 完全自包含
- 便于管理

**模式 3：错误提示（安全网）**
```bash
#!/bin/bash
osascript -e 'display alert "PyMOL Not Found" message "Please install PyMOL.app from https://pymol.org/ or run the installer again to install conda PyMOL."'
exit 1
```

**作用：**
- 防止静默失败
- 给用户明确指引

### 3. 新增文档

#### `INSTALL.md`
包含：
- 问题诊断
- 3 种解决方案
- 为什么会出问题
- 验证和启动方式
- 常见问题解答

## 技术对比

### 旧版 vs 新版

| 特性 | 旧版 | 新版 |
|------|------|------|
| PyMOL 检测 | ❌ 不检测 | ✅ 智能检测 |
| 安装验证 | ❌ 无验证 | ✅ 验证成功 |
| 系统 PyMOL | ❌ 忽略 | ✅ 优先使用 |
| 错误处理 | ❌ 静默失败 | ✅ 明确提示 |
| 启动器类型 | 1 种 | 3 种（自适应） |
| 用户反馈 | ❌ 无提示 | ✅ 显示类型 |

## 使用指南

### 对于用户

**情况 1：首次安装**
```bash
bash install_gluetk.sh
```

**情况 2：已安装但无法启动**
```bash
bash fix_pymol.sh
```

**情况 3：想要最佳体验**
1. 安装 PyMOL.app：https://pymol.org/
2. 运行：`bash fix_pymol.sh`

### 对于开发者

**重新打包 Installer**
```bash
# 1. 重建安装器 App（已修复）
bash rebuild_installer.sh

# 2. 打包 DMG
bash package_dmg.sh
```

现在生成的 Installer 会：
- ✅ 检测系统 PyMOL.app
- ✅ 验证 conda PyMOL 安装
- ✅ 生成正确的启动器
- ✅ 给出明确提示

## 验证步骤

### 验证修复是否成功

```bash
# 1. 检查脚本语法
bash -n rebuild_installer.sh
bash -n install_gluetk.sh
bash -n fix_pymol.sh

# 2. 测试快速修复
bash fix_pymol.sh

# 3. 检查生成的启动器
cat ~/Desktop/GlueTK.app/Contents/MacOS/launcher

# 4. 尝试启动
open ~/Desktop/GlueTK.app
```

### 预期结果

- **有 PyMOL.app**：启动器使用系统 PyMOL.app
- **只有 conda PyMOL**：启动器使用 conda run
- **都没有**：显示错误对话框，不会静默失败

## 影响范围

### 修改的文件
1. ✅ `rebuild_installer.sh` - 修复启动器生成逻辑
2. ✅ `install_gluetk.sh` - 新增完整安装脚本
3. ✅ `fix_pymol.sh` - 新增快速修复脚本
4. ✅ `INSTALL.md` - 新增安装文档
5. ✅ `CHANGELOG_INSTALLER_FIX.md` - 本文件

### 不影响的部分
- ❌ `gluetk/` 源码（无需修改）
- ❌ `package_dmg.sh`（打包脚本）
- ❌ PyMOL 插件逻辑

## 后续建议

### 短期
1. ✅ 测试修复脚本
2. 🔄 重新打包 DMG
3. 🔄 更新发布说明

### 长期
1. 考虑提供预编译的 PyMOL bundle
2. 添加更多错误恢复机制
3. 支持其他平台（Linux, Windows）

## 已知限制

1. **macOS 专用**：当前修复只针对 macOS
2. **PyMOL.app 路径固定**：假设安装在 `/Applications/PyMOL.app`
3. **conda 依赖**：仍然需要 conda（但不再强制 conda PyMOL）

## 测试建议

### 测试场景

**场景 1：干净系统 + PyMOL.app**
- 安装 PyMOL.app
- 运行 `install_gluetk.sh`
- 预期：使用 PyMOL.app

**场景 2：干净系统 + 无 PyMOL**
- 不安装 PyMOL.app
- 运行 `install_gluetk.sh`
- 预期：安装 conda PyMOL

**场景 3：旧安装修复**
- 已有 gluetk 环境但无 PyMOL
- 运行 `fix_pymol.sh`
- 预期：安装 PyMOL 并重建启动器

**场景 4：完全失败情况**
- conda PyMOL 安装失败
- 没有 PyMOL.app
- 预期：显示错误对话框

---

**修复完成！✨**

所有改进已实现，现在的安装器更健壮、更智能、更用户友好。
