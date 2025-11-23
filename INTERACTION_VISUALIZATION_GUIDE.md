# 相互作用可视化指南 (Interaction Visualization Guide)

## 相互作用线（虚线）的含义

### 颜色编码 (Color Coding)
每种颜色的虚线代表不同类型的分子间相互作用：

| 颜色 | 相互作用类型 | 英文 | 重要性 | 典型距离 |
|------|-------------|------|--------|---------|
| 🟡 **黄色** | 氢键 | Hydrogen Bond | ⭐⭐⭐⭐⭐ | 2.5-3.5 Å |
| 🟣 **品红色** | 盐桥 | Salt Bridge | ⭐⭐⭐⭐⭐ | 3.0-4.5 Å |
| 🟢 **绿色** | 疏水相互作用 | Hydrophobic | ⭐⭐⭐ | 3.5-4.5 Å |
| 🟠 **橙色** | π-π 堆积 | π-π Stacking | ⭐⭐⭐⭐ | 3.3-5.5 Å |
| 🔵 **青色** | 卤素键 | Halogen Bond | ⭐⭐⭐⭐ | 3.0-4.0 Å |
| 🟪 **紫色** | 金属配位 | Metal Coordination | ⭐⭐⭐⭐ | 2.0-3.0 Å |

### 相互作用的药物设计意义

1. **氢键 (黄色线)**
   - 最重要的药物-靶标相互作用
   - 提供方向特异性
   - 每个氢键贡献 ~1-2 kcal/mol 结合能

2. **盐桥 (品红色线)**
   - 带电基团间的静电吸引
   - 非常强的相互作用 (3-5 kcal/mol)
   - 对pH敏感

3. **疏水相互作用 (绿色线)**
   - 非极性基团间的接触
   - 熵驱动的相互作用
   - 贡献整体结合亲和力

4. **π-π 堆积 (橙色线)**
   - 芳香环之间的相互作用
   - 常见于药物分子与Phe/Tyr/Trp

## 如何优化显示

### 1. 只显示关键相互作用（推荐）
```python
# 只显示氢键和盐桥（最重要的相互作用）
analyze_protein_ligand_interactions('structure', 'LIG', 
                                   key_interactions_only=True)
```

### 2. 控制显示的残基
```python
# 使用新的参数只显示有相互作用的残基
highlight_csv_residues('interactions.csv', obj='structure',
                      show_only_interactions=True,  # 只显示相互作用残基
                      show_labels=True,             # 显示残基标签
                      show_interaction_type=True)   # 标签中包含相互作用类型
```

### 3. 调整虚线样式
在 PyMOL 命令行中：
```
# 让虚线更明显
set dash_width, 3      # 加粗虚线
set dash_gap, 0.3      # 减少间隙
set dash_length, 0.5   # 增加虚线段长度

# 隐藏距离数值（如果太乱）
hide labels, dist_*
hide labels, interact_*
```

### 4. 清理视图
```python
# 隐藏所有，只显示关键部分
hide everything
show cartoon, polymer.protein
show sticks, intsel_*    # 只显示相互作用残基
show dashes              # 显示相互作用线
```

## 标签含义

标签格式：`残基名+残基号`，例如：
- `R123` = 精氨酸(Arginine) 123号
- `D85` = 天冬氨酸(Aspartic acid) 85号
- `Y67` = 酪氨酸(Tyrosine) 67号

如果启用 `show_interaction_type=True`，标签会变成：
- `R123(SB)` = 精氨酸123号，参与盐桥(Salt Bridge)
- `S45(HB)` = 丝氨酸45号，参与氢键(Hydrogen Bond)

## 实用命令

### 查看特定类型的相互作用
```python
# 只显示氢键
hide dashes
show dashes, dist_*HB*

# 只显示盐桥
hide dashes  
show dashes, dist_*SB*
```

### 导出高质量图片
```python
# 设置高质量渲染
set ray_trace_mode, 1
set surface_quality, 2
bg_color white
ray 2400, 2400
png interaction_figure.png, dpi=300
```

## 常见问题

### Q: 为什么显示了很多没有标签的残基？
A: 默认会显示所有蛋白链。使用 `show_only_interactions=True` 参数只显示参与相互作用的残基。

### Q: 虚线太多看不清？
A: 可以：
1. 使用 `key_interactions_only=True` 只分析关键相互作用
2. 在 PyMOL 中手动隐藏某些类型：`hide dashes, *Hydrophobic*`

### Q: 如何知道哪条线连接哪两个原子？
A: 点击 PyMOL 中的虚线对象，会在控制台显示连接的两个原子信息。

### Q: 配体和蛋白的颜色如何区分？
- **配体**：橙色碳原子 (tv_orange)
- **A链蛋白**：浅蓝色碳原子 (lightblue)
- **B链蛋白**：浅橙色碳原子 (lightorange)
- **其他元素**：N=蓝色，O=红色，S=黄色

---

## 快速示例

完整的分析和可视化流程：

```python
# 1. 加载结构
fetch 1hsg  # 或 load your_structure.pdb

# 2. 分析相互作用（只分析关键的）
analyze_protein_ligand_interactions('1hsg', 'MK1', 
                                   key_interactions_only=True,
                                   output_csv='key_interactions.csv')

# 3. 可视化（只显示相互作用残基）
highlight_csv_residues('key_interactions.csv', obj='1hsg',
                      show_only_interactions=True,
                      show_labels=True,
                      show_atom_lines=True)

# 4. 优化显示
set cartoon_transparency, 0.7
zoom intsel_*
orient

# 5. 查看结果
# 黄色虚线 = 氢键（最重要）
# 品红虚线 = 盐桥（很重要）
# 绿色虚线 = 疏水（如果显示的话）
```

这样您就能得到一个清晰的、只显示关键相互作用的图像。