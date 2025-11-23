# GlueTK Plugin for PyMOL

蛋白质相互作用分析与可视化插件

## 🚀 快速开始

### 基本分析

```python
# 加载结构
load your_complex.pdb

# 分析相互作用（默认标准）
analyze_protein_ligand_interactions(
    obj_name='your_complex',
    ligand_resname='LIG'
)
```

### 使用 Schrödinger 标准（严格模式）

```python
# Schrödinger 标准：氢键 ≤2.8Å，更严格
analyze_protein_ligand_interactions(
    obj_name='your_complex',
    ligand_resname='LIG',
    use_schrodinger_standard=True  # 启用 Schrödinger 标准
)
```

## 📊 两种标准对比

| 相互作用 | 默认标准 | Schrödinger 标准 |
|---------|---------|-----------------|
| 氢键 | ≤ 3.5 Å | ≤ 2.8 Å (Maestro) |
| 盐桥 | ≤ 4.0 Å | ≤ 4.0 Å |
| 疏水 | ≤ 4.5 Å | ≤ 4.5 Å |
| 适用场景 | 快速筛选 | 发表级别分析 |

## 💡 主要功能

### 1. 相互作用分析

```python
# 基本分析
analyze_pdb_interactions('protein', output_csv='results.csv')

# 蛋白-配体分析
analyze_protein_ligand_interactions('protein', 'LIG')

# 三元复合体分析
analyze_ternary_complex('protein', 'PROTAC', ['A'], ['B'])

# 原子对分析
analyze_atom_pair_interactions('protein', 'elem N', 'elem O', 3.5)
```

### 2. 可视化

```python
# 3D 可视化
visualize_protein_ligand_3d('protein', result, 'LIG')

# 控制相互作用线条显示
toggle_interaction_lines(True, 3.0)   # 显示粗线条
toggle_interaction_lines(False)       # 隐藏所有线条

# 生成网络图
generate_interaction_network_plot(result, output_path='network.png')

# 高亮显示
highlight_csv_residues('interactions.csv', obj='protein')
```

### 3. GUI 界面

```python
# 打开 GUI
gluetk_gui

# 蛋白-配体专用 GUI
protein_ligand_gui
```

## 📖 详细示例

### 示例 1：HIV-1 蛋白酶分析

```python
# 获取结构
fetch 1hsg

# 默认标准分析（快速）
analyze_protein_ligand_interactions(
    obj_name='1hsg',
    ligand_resname='MK1',
    output_csv='1hsg_default.csv'
)

# Schrödinger 标准分析（严格）
analyze_protein_ligand_interactions(
    obj_name='1hsg',
    ligand_resname='MK1',
    output_csv='1hsg_schrodinger.csv',
    use_schrodinger_standard=True
)
```

### 示例 2：批量分析

```python
ligands = ['ATP', 'NAD', 'GDP']

for lig in ligands:
    analyze_protein_ligand_interactions(
        obj_name='protein',
        ligand_resname=lig,
        output_csv=f'{lig}_interactions.csv',
        use_schrodinger_standard=True  # 使用严格标准
    )
```

### 示例 3：导出高质量图片

```python
# 分析
result = analyze_protein_ligand_interactions('protein', 'LIG', 
                                             use_schrodinger_standard=True)

# 可视化
visualize_protein_ligand_3d('protein', result, 'LIG')

# 调整相互作用线条
toggle_interaction_lines(True, 3.5)  # 显示更粗的线条
set dash_color, marine, interact_Hbond_*  # 改变氢键颜色

# 设置渲染
bg_color white
set ray_opaque_background, 1
set ambient, 0.4
set spec_power, 200

# 导出
png publication.png, width=3000, height=2000, dpi=300, ray=1
```

## 🔧 参数说明

### Schrödinger 标准参数（内置）

```python
SCHRODINGER_PARAMS = {
    "hbond": {"max_distance": 2.8},        # Maestro 标准
    "saltbridge": {"max_distance": 4.0},
    "hydrophobic": {"max_distance": 4.5},
    "pi_pi": {"max_distance": 5.5, "min_distance": 3.3},
    "pi_cation": {"max_distance": 6.0},
}
```

### 自定义参数

如需修改参数，编辑 `interaction_analyzer.py` 中的 `SCHRODINGER_PARAMS` 字典。

## 📝 输出格式

### CSV 文件

```csv
Ligand_Chain,Ligand_Residue,Ligand_Atom,Protein_Chain,Protein_Residue,Protein_Atom,Distance,Interaction
L,LIG 301,N1,A,SER 99,OG,2.75,氢键
L,LIG 301,O2,A,ARG 120,NH1,3.85,盐桥
L,LIG 301,ring,A,PHE 45,ring,4.20,π–π 堆积
```

### 返回字典

```python
{
    'ligand_residues': [...],
    'protein_chains': [...],
    'interactions': [...],
    'standard_used': 'Schrödinger',  # 或 'Default'
    'parameters': {
        'hbond_cutoff': 2.8,
        'saltbridge_cutoff': 4.0,
        ...
    }
}
```

## ❓ 常见问题

### Q: 如何选择合适的标准？

**A:** 
- **默认标准**：快速筛选、初步分析
- **Schrödinger 标准**：发表论文、精确分析

### Q: 找不到配体？

**A:** 
```python
# 方法1：使用 organic 自动检测
ligand_selection='organic'

# 方法2：查看所有残基
iterate all, print(resn)

# 方法3：显示非标准残基
show sticks, hetatm and not solvent
```

### Q: 氢键数量为 0？

**A:** 
```python
# 添加氢原子
h_add

# 然后重新分析
analyze_protein_ligand_interactions(...)
```

### Q: Schrödinger 标准结果太少？

**A:** 
- 这是正常的，Schrödinger 标准更严格（氢键 ≤2.8Å vs 默认 ≤3.5Å）
- 适合发表级别的精确分析
- 如需更多结果，使用默认标准

## 📚 获取帮助

```python
# 查看函数帮助
help(analyze_protein_ligand_interactions)

# 查看所有命令
help(gluetk)

# 打开 GUI
gluetk_gui
```

## 🎯 参数对比表

### analyze_protein_ligand_interactions 参数

| 参数 | 类型 | 默认值 | 说明 |
|-----|------|--------|------|
| obj_name | str | None | PyMOL 对象名称 |
| ligand_resname | str | None | 配体残基名（None=自动检测） |
| protein_chains | list | None | 蛋白链列表（None=自动检测） |
| output_csv | str | None | CSV 输出路径 |
| distance_cutoff | float | 4.5 | 距离截断值（Å） |
| **use_schrodinger_standard** | **bool** | **False** | **启用 Schrödinger 标准** |

## 📄 引用

如在研究中使用，请引用：

```
GlueTK Plugin for PyMOL
蛋白-配体相互作用分析（支持 Schrödinger 标准）
版本: 1.0.0
```

参考 Schrödinger Maestro 标准：
- [https://www.schrodinger.com/](https://www.schrodinger.com/)

---

**祝您研究顺利！** 🎉

