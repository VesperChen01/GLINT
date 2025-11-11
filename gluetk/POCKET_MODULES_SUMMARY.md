# GlueTK 口袋分析模块功能总结

今天更新的口袋相关脚本功能整理（2025-11-11）

---

## 📦 模块概览

### 1. **pocket_detector.py** - 口袋检测核心模块
**最后更新**: 2025-11-11 10:01

#### 主要功能
- ✅ 基于网格法的蛋白口袋检测（纯 Python 实现）
- ✅ 几何性质分析：体积、表面积、深度、开口大小、球形度
- ✅ 化学性质分析：疏水性、极性、电荷、芳香性、氢键供受体
- ✅ 可成药性评分（druggability score）
- ✅ 口袋对比分析（比较两个结构的口袋变化）
- ✅ CSV 导出功能

#### 核心类与函数

**类**: `PocketDetector`
- 参数配置：
  - `grid_spacing`: 网格间距 (默认 0.6 Å)
  - `probe_radius`: 探针半径 (默认 1.4 Å，模拟水分子)
  - `min_volume`: 最小口袋体积 (默认 20.0 Ų)
  - `min_depth`: 最小埋藏深度 (默认 2.0 Å)
  - `max_solvent_access`: 最大溶剂可及度 (0-1)

**主要函数**:
```python
# 检测口袋
detect_pockets(obj_name=None, pdb_file=None, selection='all', 
               grid_spacing=0.6, min_volume=20.0, output_csv=None)
→ 返回: List[dict] 口袋列表

# 对比两个结构的口袋
compare_pockets(obj_a, obj_b, align=True, output_csv=None, **kwargs)
→ 返回: (pockets_a, pockets_b, comparison)
```

#### 口袋数据结构
每个口袋包含以下信息：
```python
{
    'id': int,                          # 口袋编号
    'volume': float,                    # 体积 (Ų)
    'surface_area': float,              # 表面积 (Ų)
    'depth': float,                     # 深度 (Å)
    'mouth_size': float,                # 开口大小 (Å)
    'sphericity': float,                # 球形度 (0-1)
    'solvent_access': float,            # 溶剂可及度 (0-1)
    'center': (x, y, z),                # 中心坐标
    'grid_coords': [(i,j,k), ...],      # 网格坐标
    'grid_points': [(x,y,z), ...],      # 实际坐标
    'residues': [                       # 口袋残基
        {'chain': str, 'resn': str, 'resi': str}, ...
    ],
    # 化学性质
    'hydrophobicity': float,            # 疏水性 (0-1)
    'polarity': float,                  # 极性 (0-1)
    'net_charge': int,                  # 净电荷
    'aromaticity': float,               # 芳香性 (0-1)
    'hbond_donors': int,                # 氢键供体数
    'hbond_acceptors': int,             # 氢键受体数
    'druggability_score': float         # 可成药性 (0-1)
}
```

---

### 2. **pocket_visualizer.py** - 口袋可视化模块
**最后更新**: 2025-11-11 10:04

#### 主要功能
- ✅ PyMOL CGO 对象生成（球体、网格、表面）
- ✅ 按性质着色（体积、疏水性、可成药性、静电势等）
- ✅ 口袋对比可视化（显示变化）
- ✅ 与 APBS 静电势叠加
- ✅ 与相互作用网络联动
- ✅ 标签显示

#### 核心函数

```python
# 主可视化函数
visualize_pockets(pockets, obj_name='pockets', 
                 color_by='volume',              # 着色依据
                 show_spheres=True,              # 显示球体
                 show_surface=False,             # 显示表面
                 show_mesh=False,                # 显示网格
                 sphere_radius=1.0,              # 球体半径
                 transparency=0.5)               # 透明度

# 显示标签
show_pocket_labels(pockets, obj_name='pocket_labels')

# 口袋对比可视化
visualize_pocket_comparison(comparison, pockets_a, pockets_b, 
                            obj_a_name='obj_a', obj_b_name='obj_b')

# 叠加静电势
overlay_pocket_electrostatics(pocket, apbs_map_file=None, obj_name='protein')

# 高亮口袋内的相互作用
highlight_pocket_interactions(pocket, interaction_csv=None)

# 导出口袋为 PDB 伪原子
export_pocket_to_pdb(pocket, output_pdb, grid_data=None)
```

#### 颜色方案
- **color_by** 选项：
  - `'volume'`: 按体积 (彩虹色谱)
  - `'hydrophobicity'`: 按疏水性 (白-橙)
  - `'druggability'`: 按可成药性 (灰-绿)
  - `'charge'`: 按电荷 (蓝-白-红)
  - `'depth'`: 按深度 (彩虹色谱)

- **对比颜色**：
  - 蓝色: 匹配的口袋
  - 绿色: 扩大的口袋 (ΔV > 10 Ų)
  - 黄色: 缩小的口袋 (ΔV < -10 Ų)
  - 橙色: 新增的口袋
  - 红色: 消失的口袋

---

### 3. **pocket_docking.py** - 口袋导向对接模块
**最后更新**: 2025-11-11 10:42

#### 主要功能
- ✅ 从口袋自动生成 Vina config.txt
- ✅ 支持用户自定义 config 文件
- ✅ 批量对接到多个口袋
- ✅ 整合相互作用分析
- ✅ 结果可视化

#### 核心函数

```python
# 计算口袋对接盒子参数
calculate_pocket_box(pocket_data: Dict, padding: float = 5.0)
→ 返回: {'center_x', 'center_y', 'center_z', 'size_x', 'size_y', 'size_z'}

# 生成 Vina 配置文件
generate_vina_config(receptor_pdbqt, ligand_pdbqt, box_params, 
                    output_pdbqt, config_path=None,
                    exhaustiveness=8, num_modes=9)
→ 返回: config_path

# 对接到单个口袋
dock_to_pocket(receptor_pdbqt, ligand_pdbqt, pocket_data, output_dir,
              pocket_id=1, padding=5.0, exhaustiveness=8,
              custom_config=None, vina_bin=None)
→ 返回: {'success', 'output_pdbqt', 'config_file', 'affinity', 'log'}

# 批量对接到多个口袋
dock_to_multiple_pockets(receptor_pdbqt, ligand_pdbqt, pockets, output_dir,
                        max_pockets=3, padding=5.0, exhaustiveness=8)
→ 返回: List[dict] 对接结果列表

# PyMOL 命令：口袋对接
pocket_based_docking(obj_name, ligand_file, output_dir=None,
                    auto_detect_pockets=True, custom_config=None,
                    max_pockets=3, exhaustiveness=8)
→ 返回: {'success', 'results', 'output_dir', 'n_pockets'}

# 可视化对接结果
visualize_docking_result(obj_name, docking_pdbqt, pocket_id=1, 
                        show_interactions=True)
```

#### 使用示例

```python
# 在 PyMOL 中
# 1. 自动检测口袋并对接
pocket_based_docking('protein', 'ligand.mol2', max_pockets=3)

# 2. 使用自定义配置对接
pocket_based_docking('protein', 'ligand.mol2', 
                    custom_config='my_config.txt',
                    auto_detect_pockets=False)

# 3. 可视化结果
visualize_docking_result('protein', 'output/ligand_pocket1_out.pdbqt', 
                        pocket_id=1, show_interactions=True)
```

---

### 4. **pocket_glue_integration.py** - 口袋与分子胶联动模块
**最后更新**: 2025-11-11 10:18

#### 主要功能
- ✅ 口袋与 PPI 界面关联分析
- ✅ 口袋与相互作用联动
- ✅ 三元复合体口袋对比（有无分子胶）
- ✅ 与 APBS 静电势叠加分析
- ✅ G-motif 与口袋综合分析
- ✅ 综合评分与可视化

#### 核心函数

```python
# 分析 PPI 界面上的口袋
analyze_pockets_in_ppi_interface(obj_name, chain_a, chain_b,
                                output_csv=None, visualize=True,
                                **pocket_kwargs)
→ 返回: List[dict] 界面口袋列表（添加了 interface_overlap 字段）

# 分析分子胶对口袋的影响
analyze_pockets_with_glue(obj_with_glue, obj_without_glue=None,
                         glue_chain=None, 
                         protein_a_chain=None, protein_b_chain=None,
                         output_dir=None, visualize=True,
                         **pocket_kwargs)
→ 返回: (pockets_with, pockets_without, comparison)

# 口袋与相互作用关联
correlate_pockets_with_interactions(pockets, interaction_csv, output_csv=None)
→ 返回: List[dict] 关联结果
  # 每个结果包含:
  - pocket_id: 口袋编号
  - num_interactions: 相互作用数量
  - interaction_types: {'H-bond': 5, 'Hydrophobic': 3, ...}
  - interaction_density: 相互作用密度 (/Ų)

# 口袋与 APBS 静电势整合
integrate_pockets_with_electrostatics(pockets, obj_name, 
                                     apbs_dx_file=None,
                                     output_csv=None, visualize=True)
→ 返回: List[dict] 静电势分析结果

# 一键式综合分析
comprehensive_glue_pocket_analysis(obj_with_glue, obj_without_glue=None,
                                  protein_a_chain=None, protein_b_chain=None,
                                  interaction_csv=None, apbs_dx_file=None,
                                  output_dir='glue_pocket_analysis',
                                  **pocket_kwargs)
→ 返回: 综合分析结果字典

# G-motif 与口袋综合分析（针对 CRBN）
comprehensive_gmotif_pocket_analysis(obj_name, e3_chain, substrate_chain,
                                    glue_chain=None,
                                    output_dir='gmotif_pocket_analysis')
→ 返回: 综合分析结果字典
```

#### 综合分析工作流

**综合分子胶口袋分析** (`comprehensive_glue_pocket_analysis`):
1. 检测并对比口袋（有无分子胶）
2. 关联口袋与相互作用
3. 整合静电势分析
4. 生成综合报告和可视化

**G-motif 口袋综合分析** (`comprehensive_gmotif_pocket_analysis`):
1. 识别 G-motif
2. 检测 G-motif 周围口袋
3. 分析分子胶与 G-motif 口袋的相互作用
4. 关联口袋与相互作用
5. 分析 PPI 界面口袋
6. 可视化所有结果
7. 生成综合报告

---

## 🎯 GUI 集成建议

### 当前已集成的功能（unified_gui.py）
- ✅ 口袋检测基础功能
- ✅ 参数设置（grid spacing, min volume）
- ✅ 颜色选择（Volume / Druggability / Hydrophobicity / Depth）
- ✅ Vina 对接（支持自定义 config）
- ✅ 快速评分

### 建议新增的 GUI 功能

#### 1. **高级口袋分析选项卡**
```
📊 Advanced Pocket Analysis
├── 口袋对比 (Pocket Comparison)
│   ├── Object A: [dropdown]
│   ├── Object B: [dropdown]
│   ├── Align structures: [checkbox]
│   └── [Compare Pockets] [Visualize Comparison]
│
├── 界面口袋分析 (PPI Interface Pockets)
│   ├── Chain A: [input]
│   ├── Chain B: [input]
│   └── [Analyze Interface Pockets]
│
└── G-motif 口袋分析 (G-motif Pocket Analysis)
    ├── E3 Chain: [input]
    ├── Substrate Chain: [input]
    ├── Glue Chain: [input] (optional)
    └── [Comprehensive Analysis]
```

#### 2. **口袋-相互作用关联**
```
🔗 Pocket-Interaction Correlation
├── Detected Pockets: [table showing pockets]
├── Interaction CSV: [file picker]
└── [Correlate] [Show Results]

结果显示：
- Pocket ID
- Volume
- Num Interactions
- Interaction Density
- Druggability Score
```

#### 3. **口袋静电势分析**
```
⚡ Pocket Electrostatics
├── APBS .dx file: [file picker]
├── Selected Pocket: [dropdown]
└── [Overlay Electrostatics]
```

#### 4. **批量对接到口袋**
```
🎯 Batch Pocket Docking
├── Detected Pockets: [table - checkboxes]
├── Ligand: [file picker]
├── Max Pockets: [slider 1-10, default 3]
├── Exhaustiveness: [slider 1-32, default 8]
└── [Dock to Selected Pockets]

结果显示：
- Pocket ID
- Affinity (kcal/mol)
- Output File
- [Load] button for each
```

---

## 📊 输出文件格式

### 1. 口袋检测 CSV
```csv
Pocket_ID,Volume_A3,Surface_Area_A2,Depth_A,Mouth_Size_A,Sphericity,
Solvent_Access,Center_X,Center_Y,Center_Z,Hydrophobicity,Polarity,
Net_Charge,Aromaticity,HBond_Donors,HBond_Acceptors,Druggability_Score,Residues
1,45.2,32.1,5.3,8.2,0.65,0.12,10.5,20.3,15.8,0.68,0.25,1,0.15,3,5,0.72,"A:TRP:123;A:PHE:145;..."
```

### 2. 口袋对比 CSV
```csv
Match_Type,Pocket_A_ID,Pocket_B_ID,Distance_A,Delta_Volume_A3,
Delta_Surface_Area_A2,Delta_Druggability
matched,1,1,1.2,5.3,3.1,0.05
new,,2,,25.0,18.0,0.68
lost,3,,,,-22.0,-15.0,-0.15
```

### 3. 界面口袋 CSV
```csv
Pocket_ID,Volume_A3,Druggability_Score,Interface_Overlap,
Hydrophobicity,Net_Charge,Center_X,Center_Y,Center_Z,Interface_Residues
1,45.2,0.72,0.65,0.68,1,10.5,20.3,15.8,"A:123;B:456;..."
```

### 4. 口袋-相互作用关联 CSV
```csv
Pocket_ID,Volume_A3,Druggability_Score,Num_Interactions,
Interaction_Density,Interaction_Types
1,45.2,0.72,12,0.2654,"H-bond:5;Hydrophobic:4;Pi-Pi:3"
```

---

## 🚀 使用示例

### 示例 1: 基础口袋检测
```python
# 在 PyMOL 中
from pocket_detector import detect_pockets

# 检测口袋
pockets = detect_pockets('protein', 
                        grid_spacing=0.6,
                        min_volume=20.0,
                        output_csv='pockets.csv')

# 可视化
from pocket_visualizer import visualize_pockets
visualize_pockets(pockets, color_by='druggability')
```

### 示例 2: 分子胶口袋对比
```python
from pocket_glue_integration import analyze_pockets_with_glue

# 对比分析
results = analyze_pockets_with_glue(
    obj_with_glue='ternary_complex',
    obj_without_glue='binary_complex',
    protein_a_chain='A',
    protein_b_chain='B',
    output_dir='glue_analysis',
    visualize=True
)

print(f"Pockets with glue: {len(results[0])}")
print(f"Pockets without glue: {len(results[1])}")
```

### 示例 3: G-motif 口袋综合分析
```python
from pocket_glue_integration import comprehensive_gmotif_pocket_analysis

# 一键式分析
results = comprehensive_gmotif_pocket_analysis(
    obj_name='6hn0',
    e3_chain='A',           # CRBN
    substrate_chain='B',    # GSPT1
    glue_chain='L',         # Lenalidomide
    output_dir='gmotif_analysis'
)

# 查看主要结合口袋
main_pocket = results['main_binding_pocket']
print(f"Main pocket: #{main_pocket['pocket_id']}")
print(f"  Interactions: {main_pocket['num_interactions']}")
print(f"  Druggability: {main_pocket['druggability_score']:.3f}")
```

### 示例 4: 口袋导向对接
```python
from pocket_docking import pocket_based_docking

# 自动检测口袋并对接
result = pocket_based_docking(
    obj_name='protein',
    ligand_file='compound.mol2',
    max_pockets=3,
    exhaustiveness=16
)

# 可视化最佳结果
if result['success']:
    best_result = sorted(result['results'], 
                        key=lambda r: r['affinity'])[0]
    visualize_docking_result('protein', 
                            best_result['output_pdbqt'],
                            best_result['pocket_id'])
```

---

## 🔧 依赖要求

- **必需**: NumPy, SciPy, PyMOL
- **可选**: AutoDock Vina (用于对接), APBS (用于静电势)

---

## 📝 注意事项

1. **性能优化**: 
   - 网格间距越小，检测越精确，但计算时间越长
   - 大型蛋白（>1000残基）建议使用 `grid_spacing=0.8-1.0`

2. **参数调优**:
   - `min_volume`: 太小会检测到噪声口袋，太大会遗漏小口袋
   - 推荐值：15-30 Ų

3. **内存使用**:
   - 网格大小 ≈ (蛋白尺寸 / grid_spacing)³
   - 超大蛋白可能需要较大内存

4. **对接注意**:
   - 确保 Vina 已安装: `conda install -c conda-forge vina`
   - PDBQT 格式要求正确的原子类型和电荷

---

## 📚 参考文献

口袋检测算法参考：
- FPocket: Le Guilloux, V. et al. BMC Bioinformatics (2009)
- PocketDepth: Depth-based pocket identification
- Druggability Score: Halgren, T. A. J. Chem. Inf. Model. (2009)

---

**文档版本**: 1.0  
**最后更新**: 2025-11-11  
**维护者**: Vesper Chen
