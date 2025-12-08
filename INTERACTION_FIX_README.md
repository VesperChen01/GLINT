# 🔧 相互作用检测修复说明

## 问题描述
蛋白-配体(Protein-Ligand)和蛋白-蛋白(Protein-Protein)相互作用检测都返回0个结果。

## 已应用的修复

### 1. 放宽检测参数 ✅

原始标准太严格，已放宽为更实用的参数：

| 相互作用类型 | 原始标准 | 新标准 | 说明 |
|------------|---------|--------|------|
| **氢键** | ≤3.2Å, ∠≥120° | ≤3.5Å, ∠≥100° | 更宽松的距离和角度要求 |
| **盐桥** | ≤4.5Å | ≤5.0Å | 增加距离阈值 |
| **疏水** | ≤4.0Å | ≤4.5Å | 扩大疏水接触范围 |
| **检测优先级** | 氢键 > 盐桥 > 疏水 | 盐桥 > 氢键 > 疏水 | 修正为标准顺序，保持互斥 |

**修改文件**: `gluetk/interaction_analyzer.py` (第82-95行)

### 2. 启用调试打印 ✅

现在运行分析时会显示：
- ✅ 检测到的配体和蛋白质残基数量
- 🔬 使用的检测参数（距离、角度阈值）
- ✅ 找到的相互作用数量
- ⚠️ 如果检测到0个相互作用，会给出可能原因和建议

**修改文件**: `gluetk/interaction_analyzer.py` (第1506-1517行, 1721-1728行)

### 3. 修正盐桥/氢键检测顺序 ✅

**为什么盐桥和氢键要互斥？**

盐桥的本质 = **静电相互作用 + 氢键**

例子：Arg-NH₃⁺ ··· ⁻OOC-Asp
- NH₃⁺ 中的 H 原子与 COO⁻ 中的 O 原子之间**本身就形成氢键**
- 同时存在**静电吸引**（正负电荷）

**正确做法**：
- 优先识别为**盐桥**（更特异、更强）
- 不再重复报告为氢键
- 避免相互作用数量虚高

**检测优先级**（参考 PLIP, LigPlot+, Maestro）：
1. **盐桥** → 最高优先级
2. **氢键** → 次优先级
3. **疏水** → 最低优先级

**修改文件**: `gluetk/interaction_analyzer.py` (第1551-1569行)

## 如何使用修复后的版本

### 方法1: 直接在GUI中使用（推荐）

1. **重启GlueTK** - 确保加载修改后的代码
2. 加载你的蛋白-配体复合物结构
3. 在GUI的"相互作用分析"标签中：
   - 选择对象
   - 输入配体名称（或留空自动检测）
   - 点击"Analyze"

4. **查看控制台输出** - 会显示详细的检测过程：
   ```
   [analyze_protein_ligand_interactions] ✓ 检测到 1 个配体分子
   [analyze_protein_ligand_interactions] ✓ 检测到 152 个蛋白质残基 (链: A)
   [analyze_protein_ligand_interactions] 🔬 检测参数:
     • 氢键: D···A ≤ 3.5 Å, ∠D–H···A ≥ 100°
     • 盐桥: ≤ 5.0 Å
     • 疏水: ≤ 4.5 Å
     • 距离截断: 4.5 Å
   [analyze_protein_ligand_interactions] ✅ 分析完成: 找到 12 个相互作用
   ```

### 方法2: 在PyMOL命令行使用

```python
# 加载结构
load your_complex.pdb

# 运行分析（会使用新的放宽参数）
analyze_protein_ligand_interactions complex, LIG, output_csv=interactions.csv
```

### 方法3: 使用诊断工具

我创建了一个诊断脚本帮助排查问题：

```python
# 在PyMOL中运行
run /Volumes/data/git/GlueTK/diagnose_interaction.py
diagnose_interaction
```

这会显示：
- 已加载的对象
- 残基类型统计
- 检测到的配体
- 不同距离下的接触数
- 依赖检查（RDKit/NumPy等）
- 实际运行测试

## 如果仍然检测不到相互作用

如果应用修复后仍然是0个相互作用，可能原因：

### 1. 配体距离蛋白质太远
**检查方法**:
```python
# 在PyMOL中运行
distance check_dist, (resn LIG), (polymer), cutoff=4.5
```
如果没有显示任何线条，说明配体和蛋白质距离>4.5Å

**解决方案**: 增加GUI中的"Distance Cutoff"参数到5.0或6.0

### 2. 配体名称错误
**检查方法**:
```python
# 查看所有残基名称
stored.resnames = []
iterate all, stored.resnames.append(resn)
print(set(stored.resnames))
```

**解决方案**: 使用正确的配体残基名称（大小写敏感！）

### 3. 缺少RDKit依赖
**检查方法**:
```bash
python3 -c "import rdkit; print('RDKit OK')"
```

**解决方案**:
```bash
pip install rdkit numpy scipy
```

### 4. 结构文件问题
- 确保PDB文件包含ATOM/HETATM记录
- 确保链ID和残基编号正确
- 确保配体和蛋白质在同一个对象中

## 进一步放宽参数（如果需要）

如果你想要更宽松的检测标准，可以修改 `gluetk/interaction_analyzer.py` 第82-95行：

```python
INTERACTION_PARAMS = {
    "hbond": {
        "max_DA_dist": 4.0,      # 更宽松：原3.5Å
        "min_donor_angle": 90,   # 更宽松：原100°
        ...
    },
    "hydrophobic": {
        "other_max": 5.0         # 更宽松：原4.5Å
    },
    "ionic": {
        "max_dist": 5.5          # 更宽松：原5.0Å
    },
    ...
}
```

⚠️ **注意**: 过于宽松的标准会导致假阳性（检测到实际不存在的相互作用）

## 恢复原始严格标准

如果你想恢复原始的严格标准，修改相同位置为：

```python
INTERACTION_PARAMS = {
    "hbond": {
        "max_DA_dist": 3.2,
        "min_donor_angle": 120,
        ...
    },
    "hydrophobic": {
        "other_max": 4.0
    },
    "ionic": {
        "max_dist": 4.5,
        "exclude_if_hbond": True
    },
    ...
}
```

## 技术细节

### 修改的核心参数
位置: `gluetk/interaction_analyzer.py` Line 82-95

这些参数被以下模块共享使用：
- `analyze_protein_ligand_interactions()` - 蛋白-配体分析
- `analyze_protein_protein_interface()` - 蛋白-蛋白分析（via `ppi_analyzer.py`）
- `analyze_ternary_complex()` - 三元复合体分析

### 为什么原始标准太严格？

原始标准基于高分辨率晶体结构的理想几何：
- 氢键 3.2Å + 120°角度 ≈ Schrödinger Maestro标准
- 适用于≤1.5Å分辨率的结构

但实际情况：
- 大多数结构分辨率 2-3Å
- 侧链构象有不确定性
- 氢原子通常是推断的
- 结构可能有小的几何误差

**新的放宽标准更接近文献中常用的实用标准。**

## 需要帮助？

如果问题仍未解决：

1. 运行诊断脚本（上面提到的 `diagnose_interaction.py`）
2. 复制控制台的完整输出
3. 提供你的结构文件（PDB）
4. 提供配体名称和链信息

---

**修改时间**: 2025-12-08  
**修改作者**: GlueTK Maintenance  
**影响范围**: 所有相互作用检测功能
