# 2D 蛋白-配体相互作用图 - Schrödinger 风格重构

## 🎨 重构概要

已成功将 2D 相互作用图从 Discovery Studio 风格升级为 **Schrödinger Maestro** 风格，实现更简洁、现代、优雅的可视化效果。

---

## ✨ 核心改进

### 1. **配体表示优化**
- ✅ 已使用 `removeHs=True` 隐藏非极性氢原子
- ✅ 增大内边距 (`padding=0.2`)，配体更居中
- ✅ 加粗化学键 (`bondLineWidth=4`)，视觉更清晰
- ✅ 固定键长 (`fixedBondLength=30`)，结构更规整

### 2. **残基气泡重构**
- ✅ 缩小气泡半径：`45` → `28`（更紧凑）
- ✅ 缩短距离：`label_dist = 140` → `100`（更靠近配体）
- ✅ **填充颜色**：气泡内部用 25% 透明度的相互作用颜色填充
- ✅ 简化标签：`Asp 127` 格式（三字母代码 + 编号）

### 3. **环形布局算法**
- ✅ 使用极坐标均匀分布，避免斥力迭代混乱
- ✅ 角度智能调整，减少重叠
- ✅ 减少迭代次数：`5` → `3`，提升性能

### 4. **相互作用线条美化**
- ✅ **氢键 & 盐桥**：带箭头实线（显示方向性）
- ✅ **疏水相互作用**：灰色扇形弧（区域表示）
- ✅ **范德华力**：默认隐藏（可选参数开启）
- ✅ 移除距离标注（默认，可通过参数启用）

### 5. **Schrödinger 配色方案**
| 相互作用类型 | 颜色 | 样式 | 箭头 |
|-------------|------|------|------|
| 氢键 (Hydrogen Bond) | `#E91E63` 粉色 | 实线 | ✅ |
| 盐桥 (Salt Bridge) | `#2196F3` 蓝色 | 实线 | ✅ |
| π-π 堆积 (Pi-Pi Stacking) | `#4CAF50` 绿色 | 实线 | ❌ |
| π 相互作用 (Pi-Interaction) | `#FF5722` 橙红 | 实线 | ❌ |
| 疏水 (Hydrophobic) | `#9E9E9E` 灰色 | 弧形 | ❌ |
| 卤素键 (Halogen Bond) | `#00BCD4` 青色 | 实线 | ❌ |
| 范德华 (van der Waals) | `#E0E0E0` 浅灰 | 虚线 | ❌ (默认不显示) |

---

## 🚀 使用方法

### 基础调用（Schrödinger 默认风格）
```python
from gluetk.interaction_2d_plot import generate_2d_interaction_diagram

generate_2d_interaction_diagram(
    csv_path="interactions.csv",
    ligand_resname="UNK",
    obj_name="your_structure",
    output_path="output_schrodinger_style.png"
)
```

### 高级参数控制
```python
generate_2d_interaction_diagram(
    csv_path="interactions.csv",
    ligand_resname="UNK",
    obj_name="your_structure",
    output_path="output_custom.png",
    
    # 可选参数
    show_distance=True,   # 显示距离标注（默认 False）
    show_vdw=True,        # 显示范德华相互作用（默认 False）
    width=1800,           # 图片宽度（默认 1600）
    height=1400,          # 图片高度（默认 1200）
    dpi=200               # 分辨率（默认 150）
)
```

### 在 PyMOL 中使用
```python
# 在 PyMOL 命令行中
from gluetk.interaction_2d_plot import generate_2d_interaction_diagram

generate_2d_interaction_diagram(
    csv_path="/path/to/interactions.csv",
    ligand_resname="LIG",
    obj_name="protein_complex",
    output_path="/tmp/test_2d_schrodinger.png"
)
```

---

## 📊 视觉对比

### 改进前（Discovery Studio 风格）
- ❌ 所有原子可见（包括氢原子），视觉杂乱
- ❌ 大气泡（半径 45），占用空间大
- ❌ 距离标注密集，遮挡结构
- ❌ 残基分散，布局不均匀

### 改进后（Schrödinger 风格）
- ✅ 简洁骨架式，隐藏非极性氢
- ✅ 小气泡（半径 28），紧凑优雅
- ✅ 无距离标注，视觉清爽
- ✅ 环形均匀布局，更对称美观
- ✅ 氢键带箭头，方向性明确
- ✅ 疏水相互作用用灰色弧表示

---

## 🔧 技术细节

### 1. 残基气泡颜色逻辑
```python
# 气泡填充使用 25% 透明度的相互作用颜色
fill_color = mcolors.to_rgba(best_style["color"], alpha=0.25)
```

### 2. 箭头实现（氢键 & 盐桥）
```python
if style.get("arrow", False):
    arrow = FancyArrowPatch(
        (tx, ty), (ex, ey),
        arrowstyle='->,head_width=0.4,head_length=0.6',
        color=style["color"],
        linewidth=2
    )
```

### 3. 疏水弧形区域
```python
if style.get("style") == "arc":
    wedge = patches.Wedge(
        (tx, ty), ld*0.6, 
        center_angle - 25/2, 
        center_angle + 25/2,
        facecolor='#9E9E9E', 
        alpha=0.15
    )
```

---

## ✅ 验证清单

测试时请检查：
- [ ] 配体结构简洁（无非极性氢）
- [ ] 残基气泡更小、更紧凑
- [ ] 无距离标注干扰（除非手动启用）
- [ ] 疏水区用灰色弧表示
- [ ] 氢键和盐桥有方向箭头
- [ ] 整体布局对称美观
- [ ] 图例只显示实际出现的相互作用类型

---

## 📌 兼容性说明

### 向后兼容
- ✅ 旧的 `DS_STYLE` 已映射到 `SCHRODINGER_STYLE`
- ✅ 所有原有 API 调用保持不变
- ✅ 默认参数提供最佳 Schrödinger 风格

### 可选回退
如果需要显示更多细节（如距离标注），可以通过参数控制：
```python
generate_2d_interaction_diagram(
    ...,
    show_distance=True,  # 恢复距离标注
    show_vdw=True        # 显示范德华相互作用
)
```

---

## 🎯 下一步建议

### 可选增强功能
1. **溶剂暴露区域标识**：添加灰色阴影表示溶剂暴露面
2. **"吉他拨片"形状残基**：实现 Schrödinger 经典的椭圆拨片形状
3. **交互式 SVG 输出**：支持鼠标悬停显示详细信息
4. **动态布局优化**：使用力导向算法进一步优化布局

### 性能优化
- ✅ 减少迭代次数（3次斥力迭代）
- ✅ 简化图例（只显示实际类型）
- ⏸️ 可选：添加缓存机制避免重复计算

---

## 📝 更新日志

**版本 2.0 - Schrödinger 风格重构** (2024-12)
- 重构配色方案为 Schrödinger 风格
- 优化残基气泡布局（极坐标环形）
- 添加氢键/盐桥箭头支持
- 实现疏水相互作用弧形表示
- 移除默认距离标注
- 缩小气泡尺寸并紧凑化布局
- 优化配体绘图参数

**版本 1.0 - Discovery Studio 风格** (原版)
- 基础相互作用图生成
- 双分子策略（全原子映射 + 重原子绘图）
- 斥力迭代布局算法

---

## 🙏 致谢

本次重构参考了 Schrödinger Maestro 的 Ligand Interaction Diagram 设计理念，实现了更现代化的 2D 可视化风格。

---

**文档生成时间**: 2024-12-14  
**重构作者**: AI 编程伙伴 🤖
