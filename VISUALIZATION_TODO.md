# GlueTK Interaction Visualization Improvements - TODO

## Summary
需要将PPI分析中的改进应用到蛋白-配体和小分子-小分子相互作用可视化。

## 已完成 ✅
### `ppi_analyzer.py` (蛋白-蛋白相互作用)
1. ✅ 使用实际相互作用原子而非CA原子
2. ✅ 氢键、疏水接触:使用具体原子(如OG1-O, CG2-CB)
3. ✅ 盐桥:使用带电原子(如NZ-OE1)
4. ✅ π-π堆积:创建环中心伪原子
5. ✅ 阳离子-π:创建阳离子中心和环中心伪原子
6. ✅ 自定义颜色(蓝/红/绿/黄/紫)
7. ✅ 所有对象名使用英文(避免PyMOL中文问题)

## 待修复 ⏳

### 1. `highlight_residues.py` (蛋白-配体相互作用)
**文件位置**: `/Users/vesper/Desktop/git/GlueTK/GlueTK/gluetk/highlight_residues.py`

**当前问题**:
- 第351行已经使用具体原子,但颜色设置可能需要改进
- 需要确认distance对象使用`cmd.color()`和`cmd.set("dash_color")`

**需要修改**:
```python
# 第351-358行附近
cmd.distance(dist_name, sel1, sel2, cutoff=cutoff)

# 需要添加:
cmd.color(color, dist_name)  # 设置整体颜色
cmd.set("dash_color", color, dist_name)  # 设置虚线颜色
cmd.set("dash_width", 2.5, dist_name)
cmd.set("dash_gap", 0.3, dist_name)
```

**建议的颜色方案**:
```python
color_map = {
    "氢键": [64/255.0, 124/255.0, 174/255.0],  # 蓝色 (与PPI一致)
    "盐桥": [215/255.0, 92/255.0, 93/255.0],   # 红色 (与PPI一致)
    "疏水相互作用": [142/255.0, 186/255.0, 141/255.0],  # 绿色 (与PPI一致)
    "π–π 堆积": "yellow",
    "π–阳离子相互作用": "purple",
    "卤素键": "cyan",
    "水桥": "lightblue",
    "金属配位": "magenta",
}
```

### 2. `ligand_ligand_analyzer.py` (小分子-小分子相互作用)
**文件位置**: `/Users/vesper/Desktop/git/GlueTK/GlueTK/gluetk/ligand_ligand_analyzer.py`

**当前问题**:
- 第243行使用选择表达式,已经是具体原子
- 第230行使用伪原子(π-π),正确
- 但颜色设置只用了`cmd.color()`,需要添加`cmd.set("dash_color")`

**需要修改**:
```python
# 第243-255行附近
cmd.distance(name, sel1, sel2)

# 修改为:
cmd.distance(name, sel1, sel2)
col = color_map.get(inter['Type'], "white")
cmd.color(col, name)
cmd.set("dash_color", col, name)  # 添加这行
cmd.set("dash_width", 2.5, name)
cmd.set("dash_gap", 0.3, name)
```

## 实施优先级

### 高优先级 (立即修复)
1. `highlight_residues.py` - 蛋白-配体是最常用的功能

### 中优先级
2. `ligand_ligand_analyzer.py` - 小分子-小分子相对少用

## 测试计划

### 测试1: 蛋白-配体
```python
# 在PyMOL中
fetch 1a4k
highlight_csv_residues("protein_ligand_interactions.csv", show_atom_lines=True)
# 检查: 虚线颜色是否正确,距离是否<5Å
```

### 测试2: 小分子-小分子
```python
# 在PyMOL中
analyze_ligand_ligand_interactions("complex", "resn LIG", "resn DRG", visualize=True)
# 检查: 虚线颜色是否正确
```

## 注意事项

1. **PyMOL对象命名**: 确保所有对象名使用英文,避免中文字符
2. **颜色一致性**: 保持与PPI分析相同的颜色方案
3. **距离准确性**: 确保使用实际相互作用原子,距离应该<5Å(除了π-π等)
4. **向后兼容**: 修改时保持API兼容性

## 完成标准

- [ ] 所有相互作用类型使用一致的颜色
- [ ] 距离标签显示真实的原子间距离
- [ ] 虚线颜色正确显示(不是默认黄色)
- [ ] 代码中的中文注释可以保留,但对象名必须是英文
- [ ] 通过测试计划中的所有测试
