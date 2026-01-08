# GLINT 安装和使用指南

## 📝 项目重命名说明

项目已从 **GlueTK** 重命名为 **GLINT** (Glue Interface Analyzer)

- 旧名称: GlueTK
- 新名称: GLINT
- 版本: v0.1.26-beta

## ✅ 安装完成

GLINT 已成功安装到你的系统！

### 安装位置

- **Conda 环境**: `glint` (Python 3.10)
- **插件目录**: `~/.pymol/startup/glint/`
- **桌面应用**: `~/Desktop/GLINT.app`
- **PyMOL 启动脚本**: `~/.pymol/startup/01_glint.py`

## 🚀 启动方式

### 方式 1: 桌面应用（推荐）

双击桌面上的 **GLINT.app** 图标

### 方式 2: 在 PyMOL 中启动

1. 打开 PyMOL
2. 在命令行输入：
   ```python
   import glint
   glint.glint_gui()
   ```

### 方式 3: 从菜单启动

在 PyMOL 菜单中选择：**Plugins → GLINT - Molecular Glue Analyzer**

## 🔧 环境激活

如果需要在命令行中使用 GLINT 环境：

```bash
conda activate glint
```

## 📦 已安装的依赖

核心依赖：
- ✅ PyMOL (系统 PyMOL.app)
- ✅ RDKit 2025.03.6
- ✅ NumPy 1.26.4
- ✅ Pandas 2.3.3
- ✅ Matplotlib 3.10.8
- ✅ SciPy 1.15.2
- ✅ Seaborn 0.13.2
- ✅ PyQt5 5.15.11
- ✅ OpenBabel 3.1.1
- ✅ AutoDock Vina 1.2.6
- ✅ PDB2PQR 3.6.1
- ✅ HADDOCK3 (最新版)

可选依赖：
- ⚠️ Open3D (表面分析)
- ⚠️ scikit-image (图像处理)

## 🧪 验证安装

在终端中运行：

```bash
conda activate glint
python -c "import glint; print(f'GLINT version: {glint.__version__}')"
```

应该输出：
```
GLINT version: v0.1.26-beta
✅ GLINT module loaded successfully!
```

## 🔄 从旧版本迁移

如果你之前安装了 GlueTK，需要注意以下变化：

### 代码变化

| 旧代码 (GlueTK) | 新代码 (GLINT) |
|----------------|----------------|
| `import gluetk` | `import glint` |
| `gluetk.gluetk_gui()` | `glint.glint_gui()` |
| `gluetk_gui` | `glint_gui` |
| `~/.pymol/startup/gluetk/` | `~/.pymol/startup/glint/` |

### 环境变化

| 旧环境 | 新环境 |
|--------|--------|
| `conda activate gluetk` | `conda activate glint` |
| `GlueTK.app` | `GLINT.app` |

## 📚 主要功能

GLINT 提供以下分析功能：

1. **蛋白-蛋白相互作用 (PPI) 分析**
   - `ppi_analyze` - 分析蛋白-蛋白界面
   - `neo_epitope_find` - 识别新表位

2. **分子胶特异性分析**
   - `degron_annotate` - G-motif 分析
   - `find_crbn_g_motif` - CRBN degron 识别

3. **口袋检测与分析**
   - `detect_pockets` - 检测结合口袋
   - `visualize_pockets` - 可视化口袋

4. **相互作用分析**
   - `analyze_protein_ligand_interactions` - 蛋白-配体相互作用
   - `analyze_ternary_complex` - 三元复合物分析

5. **突变效应分析**
   - `analyze_mutation_effects` - 突变效应
   - `ddg_heatmap` - ΔΔG 热图

6. **批量分析**
   - `batch_gmotif` - 批量 G-motif 分析
   - `batch_ppi` - 批量 PPI 分析

## 🐛 故障排除

### 问题 1: PyMOL 中找不到 glint 模块

**解决方案**：
```python
import sys
sys.path.insert(0, '/Users/vesper/.pymol/startup')
import glint
```

### 问题 2: GUI 无法启动

**解决方案**：
1. 确认 PyQt5 已安装：
   ```bash
   conda activate glint
   python -c "import PyQt5; print('PyQt5 OK')"
   ```

2. 如果失败，重新安装：
   ```bash
   conda install -c conda-forge pyqt --force-reinstall
   ```

### 问题 3: 桌面应用无法启动

**解决方案**：
查看日志文件：
```bash
cat ~/.glint_launch.log
```

## 📖 更多信息

- **GitHub**: https://github.com/VesperChen01/GLINT
- **文档**: 查看项目目录中的 Markdown 文件
- **示例**: `examples/molecular_glue_fingerprint_demo.py`

## 🔄 更新 GLINT

要更新到最新版本：

```bash
cd /Users/vesper/Desktop/git/GlueTK
git pull origin main
bash install_glint.sh
```

---

**安装日期**: 2026-01-08  
**版本**: v0.1.26-beta  
**Python**: 3.10.15  
**PyMOL**: 系统 PyMOL.app

