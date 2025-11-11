# GlueTK GUI - 高级口袋分析功能使用指南

## 📍 位置
**Docking & Scoring** 页面 → **📊 Advanced Pocket Analysis** 卡片

该卡片包含 4 个标签页，提供完整的口袋分析工作流。

---

## 🔄 Tab 1: Pocket Comparison（口袋对比）

### 功能说明
对比两个结构的口袋变化，适用于：
- 分子胶前后对比（有无配体）
- 突变体 vs 野生型
- 不同构象状态

### 使用步骤
1. **选择对象**：
   - Object A: 选择第一个结构（如 apo 结构）
   - Object B: 选择第二个结构（如 holo 结构）
   - 点击 `R` 刷新对象列表

2. **对齐选项**：
   - ✅ **Align structures before comparison**: 建议勾选，自动对齐结构

3. **运行分析**：
   - 点击 `Compare Pockets` 开始对比
   - 等待分析完成（通常需要 30-60 秒）

4. **查看结果**：
   ```
   ✅ Comparison complete:
      obj_a: 5 pockets
      obj_b: 6 pockets
      Matched: 4
      New in obj_b: 2
      Lost from obj_a: 1
   
   📈 Top volume changes:
      Pocket 1 → 1: +15.3 Ų
      Pocket 2 → 2: -8.7 Ų
   ```

5. **可视化**：
   - 点击 `Visualize` 在 PyMOL 中显示对比结果
   - **颜色含义**：
     - 🔵 蓝色：匹配的口袋（体积变化 <10 Ų）
     - 🟢 绿色：扩大的口袋（体积增加 >10 Ų）
     - 🟡 黄色：缩小的口袋（体积减少 >10 Ų）
     - 🟠 橙色：新增的口袋（仅在 B 中）
     - 🔴 红色：消失的口袋（仅在 A 中）

### 典型应用场景
```python
# 场景 1: 分析分子胶诱导的口袋变化
# Object A = 二元复合物（E3 + 底物）
# Object B = 三元复合物（E3 + 底物 + 分子胶）
# 结果显示分子胶如何重塑界面口袋

# 场景 2: 突变体效应
# Object A = 野生型
# Object B = 突变体
# 分析突变如何影响活性位点口袋
```

---

## 🔗 Tab 2: PPI Interface（界面口袋分析）

### 功能说明
专门分析蛋白-蛋白相互作用（PPI）界面上的口袋，用于：
- 寻找分子胶结合位点
- 评估界面可成药性
- 设计 PPI 抑制剂

### 使用步骤
1. **选择对象**：
   - Object: 选择包含两个蛋白的复合物
   - 点击 `R` 刷新

2. **输入链 ID**：
   - **Chain A**: 第一个蛋白的链 ID（如 `A`）
   - **Chain B**: 第二个蛋白的链 ID（如 `B`）

3. **运行分析**：
   - 点击 `Analyze Interface Pockets`
   - 系统会：
     1. 识别 PPI 界面残基
     2. 检测所有口袋
     3. 筛选位于界面的口袋
     4. 自动可视化结果

4. **查看结果**：
   ```
   🔗 Analyzing PPI interface pockets...
      Object: ternary_complex
      Interface: Chain A + Chain B
   
   ✅ Found 3 interface pockets
   
   📊 Top interface pockets:
      Pocket 2: Vol=45.2 Ų, Overlap=75%, Drug=0.68
      Pocket 5: Vol=32.1 Ų, Overlap=60%, Drug=0.55
      Pocket 8: Vol=28.5 Ų, Overlap=45%, Drug=0.42
   ```

5. **结果解读**：
   - **Volume**: 口袋体积（越大越好容纳配体）
   - **Overlap**: 与界面的重叠度（越高越靠近界面）
   - **Drug**: 可成药性评分（>0.5 为可成药）

### 典型应用
- **分子胶筛选**：找到高 overlap + 高 drug 的口袋作为靶点
- **热点识别**：overlap >70% 的口袋通常对应热点残基
- **先导化合物设计**：根据口袋性质设计分子

---

## 🔍 Tab 3: Interactions（口袋-相互作用关联）

### 功能说明
将检测到的口袋与相互作用 CSV 文件关联，用于：
- 识别活跃的结合口袋
- 分析口袋的相互作用密度
- 优先排序潜在靶点

### 前置条件
⚠️ **必须先在 "Pocket Detection" 部分检测口袋！**

### 使用步骤
1. **检测口袋**（如果还没做）：
   - 返回上方 "🔍 Pocket Detection & Visualization"
   - 选择对象，点击 `Detect Pockets`

2. **选择相互作用文件**：
   - 点击 `Browse` 选择相互作用 CSV 文件
   - 支持的格式：
     - 蛋白-配体相互作用 CSV
     - 蛋白-蛋白相互作用 CSV
     - 任何包含残基信息的 CSV

3. **运行关联**：
   - 点击 `Correlate with Pockets`
   - 系统自动匹配相互作用到口袋

4. **查看结果**：
   ```
   🔍 Correlating pockets with interactions...
      CSV: protein_ligand_interactions.csv
      Pockets: 8
   
   ✅ Correlation complete
   
   📈 Pocket-Interaction correlation:
      Pocket 2: 12 interactions (density=0.2654/Ų)
         Types: H-bond:5, Hydrophobic:4, Pi-Pi:3
      Pocket 5: 8 interactions (density=0.1987/Ų)
         Types: H-bond:3, Hydrophobic:3, Ionic:2
      Pocket 1: 3 interactions (density=0.0889/Ų)
         Types: H-bond:2, Hydrophobic:1
   ```

5. **结果解读**：
   - **Num interactions**: 相互作用数量（越多越活跃）
   - **Density**: 相互作用密度（interactions/Ų）
     - >0.2: 高活性口袋
     - 0.1-0.2: 中等活性
     - <0.1: 低活性
   - **Types**: 相互作用类型分布

### 典型应用
```python
# 场景 1: 识别主要结合口袋
# 相互作用密度最高的口袋通常是配体主要结合位点

# 场景 2: 优化先导化合物
# 根据口袋的相互作用类型设计化合物
# 例如：H-bond多的口袋适合极性基团

# 场景 3: 变构位点发现
# 低密度但存在相互作用的口袋可能是变构位点
```

---

## ✨ Tab 4: G-motif（G-motif 口袋综合分析）

### 功能说明
针对 CRBN-G-motif 系统的一键式综合分析，包括：
1. G-motif 识别
2. G-motif 周围口袋检测
3. 分子胶-G-motif 相互作用分析
4. 口袋-相互作用关联
5. PPI 界面口袋分析
6. 综合报告生成

### 适用场景
- CRBN 分子胶研究
- G-motif 底物分析
- 三元复合物表征

### 使用步骤
1. **选择对象**：
   - Object: 包含 CRBN、底物、（可选）分子胶的复合物
   - 点击 `R` 刷新

2. **输入链 ID**：
   - **E3**: CRBN 的链 ID（通常是 `A`）
   - **Sub**: 底物的链 ID（如 GSPT1，通常是 `B`）
   - **Glue**: 分子胶的链 ID（可选，如 `L`）

3. **运行分析**：
   - 点击 `Comprehensive Analysis`
   - ⌛ **注意**：这是一个综合分析，可能需要 2-5 分钟

4. **查看结果**：
   ```
   ✨ Starting comprehensive G-motif pocket analysis...
      Object: 6hn0
      E3 chain: A
      Substrate chain: B
      Glue chain: L
   
   ⌛ This may take a few minutes...
   
   ✅ Comprehensive analysis complete!
      Output directory: gmotif_pocket_analysis/
   
   🧲 G-motif detected:
      Position: 58-65
      RMSD: 1.85 Å
   
   🔍 Pockets around G-motif: 4
   
   🎯 Main binding pocket:
      Pocket ID: 2
      Interactions: 15
      Density: 0.3125/Ų
      Druggability: 0.72
   
   🔗 G-motif interface pockets: 2
   
   📄 Check the output directory for detailed reports and CSVs
   ```

5. **输出文件**（在 `gmotif_pocket_analysis/` 目录）：
   - `gmotif_detection.csv` - G-motif 检测结果
   - `gmotif_pockets.csv` - G-motif 周围口袋
   - `gmotif_glue_interactions.csv` - 分子胶相互作用
   - `gmotif_pocket_correlations.csv` - 口袋-相互作用关联
   - `interface_pockets.csv` - 界面口袋
   - `GMOTIF_POCKET_REPORT.txt` - 综合报告

### 结果解读

#### Main Binding Pocket（主要结合口袋）
- 这是分子胶最可能结合的口袋
- 高相互作用数 + 高密度 = 强结合
- 高可成药性 = 适合小分子设计

#### G-motif Interface Pockets
- 与 G-motif 区域重叠的界面口袋
- 这些口袋对 CRBN-底物识别至关重要
- 可用于设计选择性分子胶

### 典型工作流
```python
# 完整的 CRBN 分子胶分析流程：

# 1. 加载三元复合物结构
fetch 6hn0

# 2. 在 GUI 中运行 G-motif 口袋分析
# Object: 6hn0
# E3: A (CRBN)
# Sub: B (GSPT1)
# Glue: L (Lenalidomide)

# 3. 结果应用：
# - 识别 G-motif 位置
# - 找到主要结合口袋
# - 分析分子胶如何稳定界面
# - 设计新的分子胶变体

# 4. 进一步分析：
# - 使用 Pocket Comparison 对比有无分子胶
# - 使用 Interface Pockets 分析界面变化
# - 使用 Interactions 优化相互作用
```

---

## 🎯 完整工作流示例

### 案例：分析分子胶如何重塑 PPI 界面

#### 步骤 1: 准备结构
```python
# 在 PyMOL 中加载
fetch 6bn7  # 二元复合物 (CRBN + GSPT1, no glue)
fetch 6hn0  # 三元复合物 (CRBN + GSPT1 + Lenalidomide)
```

#### 步骤 2: GUI 操作序列

**2.1 检测基础口袋**（Pocket Detection 卡片）
- Object: `6bn7`
- Grid: `0.6` Å
- Min Vol: `20` Ų
- 点击 `Detect Pockets` → 可视化

**2.2 对比分析**（Advanced → Comparison 标签）
- Object A: `6bn7` (无分子胶)
- Object B: `6hn0` (有分子胶)
- ✅ Align structures
- 点击 `Compare Pockets`
- 点击 `Visualize` 查看变化

**2.3 界面口袋分析**（Advanced → PPI Interface 标签）
- Object: `6hn0`
- Chain A: `A` (CRBN)
- Chain B: `B` (GSPT1)
- 点击 `Analyze Interface Pockets`

**2.4 G-motif 综合分析**（Advanced → G-motif 标签）
- Object: `6hn0`
- E3: `A`
- Sub: `B`
- Glue: `L`
- 点击 `Comprehensive Analysis`
- 等待完成，查看 `gmotif_pocket_analysis/` 目录

#### 步骤 3: 结果整合

从各个分析中，你将得到：
1. **Comparison**: 分子胶诱导了哪些口袋变化
2. **Interface**: 界面口袋的位置和性质
3. **G-motif**: G-motif 与结合口袋的关系
4. **总体图景**: 分子胶如何通过重塑界面口袋稳定三元复合物

---

## 💡 使用技巧

### 性能优化
- **大蛋白**（>1000 残基）：使用 `grid_spacing=0.8-1.0`
- **精细分析**：使用 `grid_spacing=0.5-0.6`
- **快速筛选**：使用 `min_volume=30-50`

### 参数调优
- **min_volume**:
  - 15-20 Ų: 捕获小口袋（如小分子结合位点）
  - 25-30 Ų: 标准口袋
  - 40-50 Ų: 仅大口袋（如肽段结合位点）

### 常见问题

**Q: "No pockets detected"**
A: 尝试：
- 减小 `min_volume`
- 增大 `grid_spacing`（更快但精度降低）
- 检查结构是否完整

**Q: "Please detect pockets first"（Interactions 标签）**
A: 先在 Pocket Detection 卡片中运行口袋检测

**Q: 分析太慢**
A: 
- 使用更大的 `grid_spacing`
- 增加 `min_volume` 过滤小口袋
- 仅分析感兴趣的区域（使用 selection）

**Q: 内存不足**
A: 
- 减小蛋白尺寸（仅保留相关域）
- 增大 `grid_spacing`
- 关闭其他占用内存的程序

---

## 📚 输出文件说明

### CSV 文件格式

#### Pocket Detection CSV
```csv
Pocket_ID,Volume_A3,Surface_Area_A2,Depth_A,Druggability_Score,Center_X,Center_Y,Center_Z,Residues
1,45.2,32.1,5.3,0.72,10.5,20.3,15.8,"A:TRP:123;A:PHE:145;..."
```

#### Pocket Comparison CSV
```csv
Match_Type,Pocket_A_ID,Pocket_B_ID,Distance_A,Delta_Volume_A3,Delta_Druggability
matched,1,1,1.2,5.3,0.05
new,,2,,25.0,0.68
lost,3,,,,-22.0,-0.15
```

#### Interface Pockets CSV
```csv
Pocket_ID,Volume_A3,Druggability_Score,Interface_Overlap,Center_X,Center_Y,Center_Z
1,45.2,0.72,0.65,10.5,20.3,15.8
```

#### Pocket-Interaction Correlation CSV
```csv
Pocket_ID,Volume_A3,Druggability_Score,Num_Interactions,Interaction_Density,Interaction_Types
1,45.2,0.72,12,0.2654,"H-bond:5;Hydrophobic:4;Pi-Pi:3"
```

---

## 🔬 高级应用

### 1. 批量分析多个结构
```python
# 使用 PyMOL 脚本自动化
structures = ['struct1', 'struct2', 'struct3']
for s in structures:
    detect_pockets(s, output_csv=f'{s}_pockets.csv')
```

### 2. 与 APBS 静电势整合
```python
# 先运行 APBS（在 Electrostatics 页面）
# 然后在口袋分析中会自动叠加静电势信息
```

### 3. 导出口袋为 PDB
```python
# 使用命令行
from pocket_visualizer import export_pocket_to_pdb
export_pocket_to_pdb(pocket, 'pocket_1.pdb')
```

---

**文档版本**: 1.0  
**最后更新**: 2025-11-11  
**适用 GlueTK 版本**: 1.0.0+
