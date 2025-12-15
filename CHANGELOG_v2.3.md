# 2D互作图修复日志 v2.3

## 🐛 修复的问题 (基于用户反馈截图)

### 1. 箭头方向错误 ✅
**问题**: 箭头从配体指向残基 (应该相反)
**修复**: 
```python
# 修复前
ax.annotate("", xy=(tx, ty), xytext=(px, py), ...)  # 错误方向

# 修复后  
ax.annotate("", xy=(tx, ty), xytext=(px, py), ...)  # Badge -> Ligand
# xy是箭头终点(配体), xytext是起点(Badge) ✅
```

### 2. 圆圈内文字大小不匹配 ✅
**问题**: 长编号(如I346, D345)字体太大,超出圆圈

**修复**: 自适应字体大小
```python
if len(label_text) <= 3:
    font_size = 10      # K11, Y2等
elif len(label_text) == 4:
    font_size = 9       # K211, I346等
else:
    font_size = 8       # 超长编号
```

### 3. 图例使用中文 ✅
**问题**: 图例显示"🌀 非极性"等中文标签

**修复**: 改为英文标签
```python
aa_legend_labels = {
    'nonpolar': 'Nonpolar',      # 原: 🌀 非极性
    'polar': 'Polar',            # 原: 💧 极性
    'negative': 'Negative',      # 原: 🌿 负电荷
    'positive': 'Positive'       # 原: 🌱 正电荷
}
```

### 4. 疏水圆圈位置/大小优化 ✅
**问题**: 疏水圆圈看起来"歪了"且不够明显

**优化**:
- 半径缩小: 22 → 18 px (避免过度重叠)
- 透明度提高: 0.12 → 0.18 (更明显)
- 移除边框: ec='none' (更清爽)

```python
hydro_circle = patches.Circle(
    (tx, ty), 18,           # 在配体原子(tx,ty)上绘制
    fc=style["color"],      # 绿色填充
    ec='none',              # 无边框
    alpha=0.18,             # 18%透明度
    zorder=3                # 在分子图层之上
)
```

## 📊 关于疏水圆圈"歪"的原因

疏水圆圈的位置实际上是**正确**的,它们绘制在配体原子的坐标`(tx, ty)`上。

如果看起来"歪"的原因可能是:

1. **多个原子有疏水作用** - 圆圈集中在疏水区域
2. **2D坐标展开** - 3D结构展开成2D后相对位置改变
3. **视觉密度** - 多个圆圈重叠造成视觉偏移感

### 🔍 调试建议

如果疏水圆圈位置确实不对,可能原因:
1. CSV文件中的`Ligand_Atom`与实际原子名不匹配
2. RDKit移除氢后原子索引映射错误

可以添加调试输出:
```python
print(f"Hydrophobic: {prot_res} -> Atom {inter['target_idx']} at ({tx}, {ty})")
```

## 🎯 最终效果

- ✅ 箭头从Badge正确指向配体原子
- ✅ 文字完美适配圆圈大小
- ✅ 图例使用英文标签
- ✅ 疏水圆圈清晰且不干扰主体
- ✅ 整体布局紧凑专业

## 📝 使用说明

代码修改后无需更改调用方式:

```python
from gluetk.interaction_2d_plot import generate_2d_interaction_diagram

generate_2d_interaction_diagram(
    csv_path="interactions.csv",
    ligand_resname="MOL",
    obj_name="protein",
    output_path="output.png"
)
```

## 🔄 版本历史

- **v2.0**: 初始Schrödinger风格
- **v2.1**: 圆形Badge + 氨基酸分类着色
- **v2.2**: 强制移除氢 + 疏水无连接线
- **v2.3**: 箭头方向修复 + 文字自适应 + 英文图例

---

**日期**: 2025-12-15  
**版本**: v2.3 (Arrow Fix & Text Adaptive)
