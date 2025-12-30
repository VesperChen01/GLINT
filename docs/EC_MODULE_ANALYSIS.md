# GlueTK EC 模块解读与结果解析指南

## 概述

EC（Electrostatic Complementarity，静电互补性）模块是 GlueTK 的核心分析功能之一，用于评估蛋白质-配体复合物中的静电匹配程度。该模块由 4 个主要文件组成：

| 文件 | 功能 |
|------|------|
| `ligand_ec_calculator.py` | 核心 EC 计算器 |
| `ec_advanced_patches.py` | 高级 EC 补丁（σ-hole、孤对电子、桥联水） |
| `ec_analysis_improved.py` | 改进的 EC 分析（依赖检查、fallback） |
| `ligand_ec_calculator_optimized.py` | 优化版 EC 计算器（向量化、进度显示） |

---

## 1. ligand_ec_calculator.py - 核心 EC 计算器

### 主要功能

这是 EC 分析的主入口，实现完整的 EC 计算工作流：

1. **PDB2PQR 集成** - 蛋白质结构准备，添加氢原子和电荷
2. **APBS 集成** - 计算蛋白质静电势（Poisson-Boltzmann 方程）
3. **配体表面采样** - 在配体表面生成采样点
4. **Gasteiger 电荷计算** - 计算配体的部分电荷
5. **EC 分数计算** - 评估静电互补性

### 核心类

#### DXGrid (Line 126)

读取和插值 APBS 输出的 OpenDX 格式势能网格。

```python
class DXGrid:
    """
    读取 APBS OpenDX 格式文件，包含：
    - 网格维度 (nx, ny, nz)
    - 原点坐标
    - 网格间距 (delta)
    - 每个网格点的势能值
    """
    def read(self, filename: str)  # 读取 DX 文件
    def interpolate(self, points: np.ndarray) -> np.ndarray  # 在给定点插值势能
```

#### LigandSurfaceSampler (Line 297)

在配体表面生成采样点，使用 Fibonacci 球面分布实现均匀采样。

```python
class LigandSurfaceSampler:
    """
    表面采样策略：
    1. 在每个原子周围放置球面点（vdW 半径 + 探针半径）
    2. 移除位于其他原子内部的点
    3. 可选：聚类/下采样以提高效率
    """
    def sample_molecule(self, mol, conformer_id=0) -> Tuple[np.ndarray, np.ndarray]
    # 返回 (surface_points, surface_normals)
```

#### GasteigerChargeCalculator (Line 441)

计算 Gasteiger 部分电荷和 Coulomb 势能。

```python
class GasteigerChargeCalculator:
    """
    使用 RDKit 的 Gasteiger 电荷实现，计算表面点的 Coulomb 势能。
    
    增强功能 (v2.0)：
    - 可选 σ-hole 虚拟点（用于 Cl/Br/I）
    - 可选孤对电子虚拟点（用于羰基 O）
    """
    def calculate_charges(self, mol) -> np.ndarray
    def calculate_potential(self, mol, points, conformer_id=0) -> np.ndarray
```

#### PDB2PQRRunner (Line 636)

PDB2PQR 接口，支持 Python API 和命令行两种方式。

```python
class PDB2PQRRunner:
    """
    PDB2PQR 添加氢原子，根据力场分配电荷和半径。
    支持 Python API 和命令行接口。
    """
    def run(self, input_pdb, output_pqr, ph=7.4, keep_chain=True, remove_water=True) -> bool
```

#### APBSRunner (Line 828)

APBS 接口，用于求解 Poisson-Boltzmann 方程计算静电势。

```python
class APBSRunner:
    """
    APBS 求解 Poisson-Boltzmann 方程计算分子周围的静电势。
    支持 Python API 和命令行接口。
    """
    def generate_input(self, pqr_file, output_prefix, ...) -> str  # 生成 APBS 输入文件
    def run(self, input_file, working_dir=None) -> Optional[str]  # 运行 APBS，返回 DX 文件路径
```

#### ECCalculator (Line 1085)

EC 分数计算核心。

```python
class ECCalculator:
    """
    EC 测量配体静电势与蛋白质势能在配体表面的匹配程度。
    
    EC_i = -2 × φ_protein × φ_ligand / (φ_protein² + φ_ligand² + ε)
    
    - 符号相反（互补）：EC > 0（蛋白提供配体表面所需）
    - 符号相同（冲突）：EC < 0（静电排斥）
    - 完美匹配：EC → +1
    """
    def calculate_ec_local(self, phi_protein, phi_ligand) -> np.ndarray
    def calculate_ec_score(self, ec_local, weights=None) -> float
    def calculate_ec_statistics(self, ec_local) -> Dict[str, float]
```

#### ECMapWriter (Line 1183)

输出 EC 可视化文件。

```python
class ECMapWriter:
    """将 EC 值写入各种格式用于可视化"""
    @staticmethod
    def write_pseudo_pdb(filename, points, values, scale=100.0, atom_name='EC')
    @staticmethod
    def write_csv(filename, points, ec_values, phi_protein=None, phi_ligand=None)
```

### 主要函数

#### calculate_ligand_ec() (Line 1258)

计算蛋白-配体复合物的 EC。

```python
def calculate_ligand_ec(
    obj_name: str = None,           # PyMOL 对象名
    ligand_resname: str = None,     # 配体残基名（如 'LIG'）
    protein_chains: List[str] = None,  # 蛋白链 ID
    output_dir: str = None,         # 输出目录
    ph: float = 7.4,                # 质子化状态的 pH
    surface_density: float = 10.0,  # 表面采样密度（点/Å²）
    visualize: bool = True,         # 是否在 PyMOL 中可视化
    pdb_file: str = None,           # 替代：使用 PDB 文件而非 PyMOL 对象
    use_sigma_holes: bool = False,  # 为 Cl/Br/I 添加 σ-hole 虚拟点
    use_lone_pairs: bool = False,   # 为羰基 O 添加孤对电子虚拟点
    keep_bridging_waters: bool = False  # 保留结构性桥联水
) -> Optional[Dict[str, Any]]
```

**返回值：**
```python
{
    'ec_score': float,           # 总体 EC 分数
    'ec_statistics': dict,       # 详细统计
    'surface_points': np.ndarray,  # 表面点坐标
    'ec_values': np.ndarray,     # 每个点的 EC 值
    'phi_protein': np.ndarray,   # 蛋白势能
    'phi_ligand': np.ndarray,    # 配体势能
    'output_files': dict,        # 生成的文件路径
    'advanced_features': dict,   # σ-hole/孤对电子使用信息
}
```

#### analyze_ternary_ec() (Line 1497)

分析三元复合物（分子胶）的 EC。

```python
def analyze_ternary_ec(
    obj_name: str = None,
    glue_resname: str = None,       # 分子胶残基名
    protein_a_chains: List[str] = None,  # 蛋白 A 链（如 E3 连接酶）
    protein_b_chains: List[str] = None,  # 蛋白 B 链（如底物）
    output_dir: str = None,
    ph: float = 7.4,
    surface_density: float = 10.0,
    visualize: bool = True,
    pdb_file: str = None
) -> Optional[Dict[str, Any]]
```

**分析的界面：**
1. EC(A-glue)：蛋白 A 势能 vs 分子胶表面
2. EC(B-glue)：蛋白 B 势能 vs 分子胶表面
3. EC(A-B)：可选的 PPI 界面分析

#### compare_ligand_ec() (Line 2094)

比较多个配体的 EC 分数（SAR 分析）。

```python
def compare_ligand_ec(
    obj_name: str,
    ligand_resnames: List[str],  # 要比较的配体残基名列表
    protein_chains: List[str] = None,
    output_dir: str = None,
    ph: float = 7.4,
    surface_density: float = 10.0
) -> Optional[Dict[str, Any]]
```

#### calculate_ec_hotspots() (Line 2358)

识别 EC 热点区域。

```python
def calculate_ec_hotspots(
    obj_name: str,
    ligand_resname: str,
    protein_chains: List[str] = None,
    output_dir: str = None,
    ph: float = 7.4,
    surface_density: float = 15.0,  # 更高密度以获得更好的热点分辨率
    hotspot_threshold: float = 0.5  # EC 阈值
) -> Optional[Dict[str, Any]]
```

#### analyze_substituent_ec_effect() (Line 2524)

分析特定取代基原子对 EC 的贡献。

```python
def analyze_substituent_ec_effect(
    obj_name: str,
    ligand_resname: str,
    substituent_atoms: List[str],  # 要分析的原子名（如 ['F1', 'F2', 'F3']）
    protein_chains: List[str] = None,
    output_dir: str = None,
    ph: float = 7.4
) -> Optional[Dict[str, Any]]
```

---

## 2. ec_advanced_patches.py - 高级 EC 补丁

### 解决的问题

标准点电荷模型无法表达某些方向性静电效应，此模块添加虚拟点来改进。

### 核心类

#### SigmaHoleGenerator (Line 107)

为卤素（Cl/Br/I）添加 σ-hole 虚拟正电点。

```python
class SigmaHoleGenerator:
    """
    σ-hole 虚拟点生成器
    
    原理：
    - 卤素原子在 C-X 键方向有电子密度缺失（σ-hole）
    - 这导致该方向呈现正电势，可与 Lewis 碱形成卤素键
    - 标准点电荷模型无法表达这种各向异性
    
    使用方法：
        generator = SigmaHoleGenerator()
        virtual_points = generator.generate(mol)
    """
    def generate(self, mol, conformer_id=0) -> Dict[str, Any]
    def generate_from_coords(self, coords, elements, bonds) -> Dict[str, Any]
```

**σ-hole 参数（基于文献值）：**

| 卤素 | 距离 (Å) | 电荷 (e) | vdW 半径 (Å) |
|------|----------|----------|--------------|
| Cl | 1.20 | +0.08 | 0.3 |
| Br | 1.30 | +0.10 | 0.3 |
| I | 1.40 | +0.12 | 0.3 |

#### LonePairGenerator (Line 299)

为羰基氧、醚氧、胺氮等添加孤对电子虚拟点。

```python
class LonePairGenerator:
    """
    孤对电子虚拟点生成器
    
    原理：
    - 羰基 C=O 的氧有两个孤对电子，位于 C=O 键两侧
    - 标准点电荷只给 O 一个负电荷，无法表达方向性
    - 添加 LP 虚拟点可改善氢键方向性预测
    """
    def generate(self, mol, conformer_id=0) -> Dict[str, Any]
```

**孤对电子参数：**

| 类型 | 距离 (Å) | 电荷 (e) | 夹角 (°) |
|------|----------|----------|----------|
| 羰基 O | 0.70 | -0.20 | 120.0 |
| 醚 O | 0.60 | -0.15 | 109.5 |
| 胺 N | 0.60 | -0.15 | 109.5 |

#### BridgingWaterFilter (Line 485)

筛选结构性重要的水分子。

```python
class BridgingWaterFilter:
    """
    桥联水筛选器
    
    从晶体结构中筛选出结构性重要的水分子：
    1. 与蛋白形成多个氢键的水
    2. 同时连接蛋白和配体的桥水
    
    原理：
    - 结构性水分子对结合亲和力有重要贡献
    - 简单删除所有水会丢失这些信息
    - 保留关键水分子可改善 EC 分析准确性
    """
    def filter_pdb(self, pdb_file, ligand_resname=None) -> Dict[str, Any]
    def write_filtered_pdb(self, input_pdb, output_pdb, waters_to_keep, ...) -> bool
```

**桥联水筛选参数：**

| 参数 | 值 | 说明 |
|------|-----|------|
| hbond_distance_max | 3.5 Å | 氢键最大距离 |
| hbond_angle_min | 120.0° | 氢键最小角度 |
| min_protein_contacts | 2 | 最少蛋白接触数 |
| ligand_distance_max | 4.0 Å | 到配体的最大距离（桥水） |

#### EnhancedChargeCalculator (Line 744)

整合所有高级补丁的增强电荷计算器。

```python
class EnhancedChargeCalculator:
    """
    增强电荷计算器
    
    整合所有高级补丁，提供统一的电荷计算接口：
    1. 基础 Gasteiger 电荷
    2. σ-hole 虚拟点
    3. 孤对电子虚拟点
    """
    def __init__(self, use_sigma_holes=True, use_lone_pairs=False, dielectric=4.0)
    def calculate(self, mol, conformer_id=0) -> Dict[str, Any]
    def calculate_potential(self, mol, points, conformer_id=0) -> np.ndarray
```

**返回值：**
```python
{
    'atom_coords': np.ndarray,      # 原子坐标 (N, 3)
    'atom_charges': np.ndarray,     # 原子电荷 (N,)
    'atom_radii': np.ndarray,       # 原子半径 (N,)
    'virtual_coords': np.ndarray,   # 虚拟点坐标 (M, 3)
    'virtual_charges': np.ndarray,  # 虚拟点电荷 (M,)
    'virtual_radii': np.ndarray,    # 虚拟点半径 (M,)
    'all_coords': np.ndarray,       # 所有点坐标 (N+M, 3)
    'all_charges': np.ndarray,      # 所有点电荷 (N+M,)
    'all_radii': np.ndarray,        # 所有点半径 (N+M,)
    'n_atoms': int,
    'n_virtual': int,
    'total_charge': float,
}
```

---

## 3. ec_analysis_improved.py - 改进的 EC 分析

### 主要改进

1. **详细的依赖检查和诊断**
2. **多层 fallback 机制**
3. **更好的错误消息和建议**
4. **性能监控和进度显示**
5. **自动降级策略**

### 核心类

#### DependencyChecker (Line 23)

依赖检查和诊断工具。

```python
class DependencyChecker:
    """检查核心依赖和 EC 专用依赖"""
    def check_core_dependencies(self) -> bool  # 检查 NumPy/SciPy/RDKit
    def check_ec_dependencies(self) -> Dict[str, bool]  # 检查 PDB2PQR/APBS
    def check_external_tools(self) -> Dict[str, bool]  # 检查命令行工具
    def get_status_report(self) -> str  # 获取状态报告
    def print_installation_guide(self)  # 打印安装指南
```

#### ECAnalysisManager (Line 147)

EC 分析管理器，处理 fallback 和错误恢复。

```python
class ECAnalysisManager:
    """管理 EC 分析环境"""
    def initialize(self) -> bool  # 初始化环境
    def get_pdb2pqr_runner(self)  # 获取 PDB2PQR 运行器
    def get_apbs_runner(self)  # 获取 APBS 运行器
    def run_ec_analysis(self, obj_name, ligand_resname, output_dir=None, **kwargs) -> Optional[Dict]
```

### Fallback 策略

| 情况 | Fallback |
|------|----------|
| PDB2PQR 不可用 | 使用简化的 PQR 转换 |
| APBS 不可用 | 使用 Coulomb 近似计算静电势 |

---

## 4. ligand_ec_calculator_optimized.py - 优化版 EC 计算器

### 性能优化

1. **向量化 Coulomb 势能计算** - 使用 NumPy 广播避免 Python 循环
2. **KDTree 加速** - 使用 SciPy cKDTree 加速距离计算
3. **分块处理** - 节省内存
4. **进度显示** - 实时显示计算进度和 ETA

### 核心类

#### ProgressTracker (Line 66)

进度跟踪和显示。

```python
class ProgressTracker:
    """跟踪和显示计算进度"""
    def __init__(self, total_steps=100, name="Processing")
    def update(self, step=None, message="")  # 更新进度
    def finish(self, message="")  # 完成进度显示
```

**显示效果：**
```
[Surface Sampling] ████████████████████░░░░░░░░░░░░░░░░░░░░  50% | 500/1000 | Elapsed: 2.5s | ETA: 2.5s
```

#### OptimizedGasteigerChargeCalculator (Line 134)

向量化 Coulomb 势能计算。

```python
class OptimizedGasteigerChargeCalculator:
    """
    优化版本的 Gasteiger 电荷计算器
    使用向量化操作加速 Coulomb 势能计算
    """
    def calculate_potential_vectorized(self, mol, points, conformer_id=0, progress_callback=None) -> np.ndarray
```

**优化策略：**
```python
# 分块处理以节省内存
chunk_size = 1000
for chunk_start in range(0, n_points, chunk_size):
    chunk_points = points[chunk_start:chunk_end]
    
    # 使用广播计算距离矩阵
    diff = chunk_points[:, np.newaxis, :] - coords[np.newaxis, :, :]
    distances = np.linalg.norm(diff, axis=2)
    
    # 向量化计算势能
    chunk_potentials = np.sum(charges[np.newaxis, :] / (dielectric * distances), axis=1)
```

#### OptimizedLigandSurfaceSampler (Line 215)

KDTree 加速的表面采样。

```python
class OptimizedLigandSurfaceSampler:
    """
    优化版本的配体表面采样器
    使用 KDTree 加速距离计算
    """
    def sample_molecule(self, mol, conformer_id=0, progress_callback=None) -> Tuple[np.ndarray, np.ndarray]
```

**优化策略：**
```python
# 构建 KDTree 以加速距离查询
if len(atoms) > 1 and SCIPY_AVAILABLE:
    atom_positions = np.array([a['pos'] for a in atoms])
    tree = cKDTree(atom_positions)
    
    # 使用 KDTree 快速过滤内部点
    distances, indices = tree.query(points, k=len(atoms))
```

#### PerformanceStats (Line 335)

性能统计。

```python
class PerformanceStats:
    """记录性能统计信息"""
    def start(self, name: str)  # 开始计时
    def end(self, name: str)  # 结束计时
    def print_summary(self)  # 打印性能摘要
```

**输出示例：**
```
============================================================
Performance Statistics
============================================================
coulomb_calculation                          5.23s ( 65.4%)
surface_sampling                             2.15s ( 26.9%)
file_io                                      0.62s (  7.7%)
------------------------------------------------------------
Total                                        8.00s (100.0%)
============================================================
```

---

## 使用示例

### 基本 EC 分析

```python
from gluetk.ligand_ec_calculator import calculate_ligand_ec

# 在 PyMOL 中
result = calculate_ligand_ec('complex', 'LIG', output_dir='./ec_output')

print(f"EC Score: {result['ec_score']:.4f}")
print(f"Positive EC fraction: {result['ec_statistics']['ec_positive_fraction']*100:.1f}%")
```

### 三元复合物（分子胶）分析

```python
from gluetk.ligand_ec_calculator import analyze_ternary_ec

result = analyze_ternary_ec(
    'complex', 
    'GLUE',           # 分子胶残基名
    ['A'],            # E3 连接酶链
    ['B'],            # 底物链
    output_dir='./ternary_ec'
)

print(f"EC(A-Glue): {result['interfaces']['A_glue']['ec_score']:.4f}")
print(f"EC(B-Glue): {result['interfaces']['B_glue']['ec_score']:.4f}")
print(f"Combined EC: {result['combined']['ec_combined_score']:.4f}")
```

### 启用 σ-hole 增强（含卤素配体）

```python
result = calculate_ligand_ec(
    'complex', 
    'LIG', 
    use_sigma_holes=True,  # 为 Cl/Br/I 添加 σ-hole
    output_dir='./ec_output'
)
```

### 多配体 SAR 比较

```python
from gluetk.ligand_ec_calculator import compare_ligand_ec

result = compare_ligand_ec(
    'complex',
    ['LIG1', 'LIG2', 'LIG3'],  # 要比较的配体
    output_dir='./ec_comparison'
)

# 输出排名
for rank, (name, score) in enumerate(result['comparison']['ranking'], 1):
    print(f"{rank}. {name}: EC = {score:.4f}")
```

### EC 热点分析

```python
from gluetk.ligand_ec_calculator import calculate_ec_hotspots

result = calculate_ec_hotspots(
    'complex',
    'LIG',
    hotspot_threshold=0.5,  # EC > 0.5 为正热点
    output_dir='./ec_hotspots'
)

print(f"Positive hotspots: {result['hotspots']['positive']['count']} points")
print(f"Negative hotspots: {result['hotspots']['negative']['count']} points")
```

### 取代基效应分析

```python
from gluetk.ligand_ec_calculator import analyze_substituent_ec_effect

result = analyze_substituent_ec_effect(
    'complex',
    'LIG',
    ['F1', 'F2', 'F3'],  # CF3 基团的三个氟原子
    output_dir='./substituent_ec'
)

print(f"Substituent EC: {result['substituent_analysis']['ec_mean']:.4f}")
print(f"Rest of ligand EC: {result['rest_of_ligand']['ec_mean']:.4f}")
```

---

## 依赖关系

### 必需依赖

| 包 | 用途 |
|-----|------|
| NumPy | 数值计算 |
| RDKit | 分子处理、Gasteiger 电荷 |

### 推荐依赖

| 包 | 用途 |
|-----|------|
| SciPy | 插值、KDTree 加速 |
| PyMOL | 可视化、结构提取 |

### EC 专用依赖

| 工具 | 用途 | 安装方式 |
|------|------|----------|
| PDB2PQR | 蛋白结构准备 | `pip install pdb2pqr` |
| APBS | 静电势计算 | `pip install apbs` 或 `brew install brewsci/bio/apbs` |

---

## 🔍 EC 结果解读详解

### 核心概念：什么是 EC？

**EC（Electrostatic Complementarity，静电互补性）** 衡量的是蛋白质和配体在结合界面的静电匹配程度。

#### 物理意义

```
蛋白质表面有正电区域 → 配体对应位置应该是负电 → 静电吸引 → EC > 0
蛋白质表面有负电区域 → 配体对应位置应该是正电 → 静电吸引 → EC > 0
蛋白质表面有正电区域 → 配体对应位置也是正电 → 静电排斥 → EC < 0
```

#### EC 计算公式

```
EC_i = -2 × φ_protein × φ_ligand / (φ_protein² + φ_ligand² + ε)
```

- `φ_protein`：蛋白质在配体表面点产生的静电势
- `φ_ligand`：配体自身在该点产生的静电势
- 当两者符号相反时，乘积为负，EC 为正（互补）
- 当两者符号相同时，乘积为正，EC 为负（冲突）

---

### 单配体 EC 分析结果解读

#### 返回的主要指标

```python
result = calculate_ligand_ec('complex', 'LIG')

# 1. 总体 EC 分数
ec_score = result['ec_score']  # 范围通常在 -1 到 +1

# 2. 详细统计
stats = result['ec_statistics']
```

#### EC 分数解读标准

| EC Score | 评价 | 说明 |
|----------|------|------|
| **> 0.5** | ⭐⭐⭐ 优秀 | 非常强的静电互补性，罕见 |
| **0.3 ~ 0.5** | ⭐⭐ 良好 | 强静电互补性，有利于结合 |
| **0.1 ~ 0.3** | ⭐ 中等 | 中等静电互补性 |
| **0 ~ 0.1** | 一般 | 弱静电互补性 |
| **< 0** | ⚠️ 不利 | 静电冲突，不利于结合 |

#### 统计指标详解

```python
stats = result['ec_statistics']

# 平均值 - 最常用的指标
ec_mean = stats['ec_mean']

# 中位数 - 对异常值更稳健
ec_median = stats['ec_median']

# 标准差 - 反映 EC 分布的均匀性
ec_std = stats['ec_std']
# 低 std (< 0.3): EC 分布均匀
# 高 std (> 0.5): EC 分布不均匀，有局部热点

# 正 EC 比例 - 最直观的指标！
positive_fraction = stats['ec_positive_fraction']
# > 60%: 大部分表面静电互补
# 50-60%: 中等
# < 50%: 静电冲突占主导

# 负 EC 比例
negative_fraction = stats['ec_negative_fraction']

# 四分位数 - 了解分布
ec_q25 = stats['ec_q25']  # 25% 分位
ec_q75 = stats['ec_q75']  # 75% 分位
```

#### 实际案例解读

**案例 1：良好的静电互补性**
```
EC Score: 0.35
EC Mean: 0.32
EC Median: 0.38
EC Std: 0.28
Positive EC fraction: 72%
Negative EC fraction: 28%

解读：
✅ EC 分数 0.35 表示良好的静电互补性
✅ 72% 的表面点显示正 EC（互补）
✅ 标准差 0.28 表示分布相对均匀
→ 结论：静电因素有利于该配体的结合
```

**案例 2：静电冲突**
```
EC Score: -0.12
EC Mean: -0.08
EC Median: -0.05
EC Std: 0.45
Positive EC fraction: 38%
Negative EC fraction: 62%

解读：
❌ EC 分数为负，表示整体静电冲突
❌ 62% 的表面点显示负 EC（冲突）
⚠️ 高标准差表示有局部差异
→ 结论：静电因素不利于结合，可能需要优化配体电荷分布
```

---

### 三元复合物（分子胶）EC 结果解读

分子胶分析会返回多个界面的 EC：

```python
result = analyze_ternary_ec('complex', 'GLUE', ['A'], ['B'])
```

#### 界面 EC 指标

```python
# 界面 1：蛋白 A（如 E3 连接酶）与分子胶
ec_a_glue = result['interfaces']['A_glue']['ec_score']

# 界面 2：蛋白 B（如底物）与分子胶
ec_b_glue = result['interfaces']['B_glue']['ec_score']

# 组合分析
combined = result['combined']
ec_combined = combined['ec_combined_score']  # 两个界面的平均
ec_asymmetry = combined['ec_asymmetry']  # 两个界面的差异
```

#### 重叠区域分析（最重要！）

```python
overlap = combined['ec_overlap']

# 重叠区域的正 EC 比例 - 分子胶效果的关键指标
pos_frac_a = overlap['pos_frac_a']  # A 界面在重叠区的正 EC 比例
pos_frac_b = overlap['pos_frac_b']  # B 界面在重叠区的正 EC 比例
pos_frac_combined = overlap['pos_frac_combined']  # 组合正 EC 比例

# 重叠区域点数
n_overlap_points = overlap['n_points']
```

#### 三元复合物解读标准

| 指标 | 优秀 | 良好 | 一般 | 不利 |
|------|------|------|------|------|
| EC(A-Glue) | > 0.3 | 0.1~0.3 | 0~0.1 | < 0 |
| EC(B-Glue) | > 0.3 | 0.1~0.3 | 0~0.1 | < 0 |
| Combined EC | > 0.3 | 0.1~0.3 | 0~0.1 | < 0 |
| Overlap Pos Frac | > 70% | 60~70% | 50~60% | < 50% |
| Asymmetry | < 0.1 | 0.1~0.2 | 0.2~0.3 | > 0.3 |

#### 三元复合物案例解读

**案例：成功的分子胶**
```
EC(A-Glue) Score: 0.28
EC(B-Glue) Score: 0.35
Combined EC: 0.315

Overlap Region (glue bridging zone - 156 points):
  Positive EC fraction (A-Glue): 68%
  Positive EC fraction (B-Glue): 75%
  Combined positive fraction: 71.5%

Asymmetry: 0.07

解读：
✅ 两个界面都有良好的 EC（0.28 和 0.35）
✅ 重叠区域 71.5% 正 EC - 分子胶桥接效果好
✅ 低不对称性（0.07）- 两个界面平衡
→ 结论：该分子胶在静电上有效地桥接两个蛋白
```

**案例：不平衡的分子胶**
```
EC(A-Glue) Score: 0.42
EC(B-Glue) Score: 0.05
Combined EC: 0.235

Overlap Region:
  Positive EC fraction (A-Glue): 78%
  Positive EC fraction (B-Glue): 52%
  Combined positive fraction: 65%

Asymmetry: 0.37

解读：
✅ A 界面 EC 很好（0.42）
⚠️ B 界面 EC 很弱（0.05）
⚠️ 高不对称性（0.37）- 两个界面不平衡
→ 结论：分子胶与蛋白 A 结合良好，但与蛋白 B 的静电互补性差
→ 建议：优化分子胶面向蛋白 B 的部分
```

---

### EC 热点分析结果解读

```python
result = calculate_ec_hotspots('complex', 'LIG', hotspot_threshold=0.5)
```

#### 热点指标

```python
hotspots = result['hotspots']

# 正热点（强互补区域）
pos_count = hotspots['positive']['count']  # 点数
pos_fraction = hotspots['positive']['fraction']  # 占比
pos_mean_ec = hotspots['positive']['mean_ec']  # 平均 EC
pos_clusters = hotspots['positive'].get('n_clusters', 0)  # 聚类数

# 负热点（冲突区域）
neg_count = hotspots['negative']['count']
neg_fraction = hotspots['negative']['fraction']
neg_mean_ec = hotspots['negative']['mean_ec']
```

#### 热点解读

| 情况 | 解读 | 建议 |
|------|------|------|
| 正热点多，负热点少 | 静电互补性好 | 保持当前设计 |
| 正热点少，负热点多 | 静电冲突严重 | 需要优化电荷分布 |
| 正负热点都集中 | 局部极化 | 检查是否有特定相互作用 |
| 正负热点都分散 | 静电效应弱 | 静电可能不是主要驱动力 |

---

### SAR 比较结果解读

```python
result = compare_ligand_ec('complex', ['LIG1', 'LIG2', 'LIG3'])
```

#### 比较指标

```python
comparison = result['comparison']

# 排名（按 EC 分数从高到低）
ranking = comparison['ranking']  # [(name, score), ...]

# 最佳/最差配体
best_ligand = comparison['best_ligand']
best_ec = comparison['best_ec_score']
worst_ligand = comparison['worst_ligand']
worst_ec = comparison['worst_ec_score']

# EC 范围
ec_range = comparison['ec_range']  # best - worst
```

#### SAR 解读

```
Rank  Ligand    EC Score    Interpretation
1     LIG2      0.38        Strong complementarity
2     LIG1      0.25        Moderate
3     LIG3      0.08        Weak

EC Range: 0.30

解读：
- LIG2 静电互补性最好，可能是活性最高的
- LIG3 静电互补性最差
- EC 范围 0.30 表示静电因素可能是 SAR 的重要驱动力
- 如果 EC 范围很小（< 0.1），静电可能不是主要因素
```

---

### 取代基效应分析结果解读

```python
result = analyze_substituent_ec_effect('complex', 'LIG', ['F1', 'F2', 'F3'])
```

#### 取代基指标

```python
sub = result['substituent_analysis']

# 取代基区域的 EC
sub_ec_mean = sub['ec_mean']
sub_pos_frac = sub['ec_positive_fraction']
sub_contribution = sub['contribution_to_total']  # 对总 EC 的贡献

# 配体其余部分的 EC
rest = result['rest_of_ligand']
rest_ec_mean = rest['ec_mean']
rest_pos_frac = rest['ec_positive_fraction']
```

#### 取代基效应解读

| 情况 | 解读 |
|------|------|
| sub_ec > rest_ec + 0.1 | 取代基改善了静电互补性 |
| sub_ec < rest_ec - 0.1 | 取代基降低了静电互补性 |
| \|sub_ec - rest_ec\| < 0.1 | 取代基对静电影响不大 |

**案例：CF3 vs CH3 效应**
```
CF3 取代基分析：
  Substituent EC Mean: 0.45
  Rest of ligand EC Mean: 0.22
  Contribution to total: 35%

解读：
✅ CF3 区域 EC（0.45）显著高于配体其余部分（0.22）
✅ CF3 贡献了 35% 的总 EC
→ 结论：CF3 的强电负性产生了有利的静电互补
→ 这可能解释为什么 CF3 类似物活性更高
```

---

### 常见问题与解答

#### Q1: EC 分数为负是否意味着配体不能结合？

**不一定。** EC 只衡量静电互补性，结合还受其他因素影响：
- 疏水相互作用
- 氢键
- 形状互补性
- 熵效应

负 EC 表示静电不利，但其他因素可能补偿。

#### Q2: 多高的 EC 才算"好"？

这取决于体系：
- 对于高度极性的结合位点，EC > 0.3 是好的
- 对于疏水性结合位点，EC 可能接近 0
- 关键是与同系列化合物比较

#### Q3: 为什么我的 EC 分数和文献不一致？

可能原因：
1. **表面采样密度不同** - 尝试调整 `surface_density` 参数
2. **介电常数不同** - 默认使用 ε=4
3. **质子化状态不同** - 检查 pH 设置
4. **是否使用 σ-hole** - 含卤素配体应启用

#### Q4: 如何改善负 EC 区域？

1. **识别负 EC 热点位置**
2. **分析该位置的蛋白静电势**
   - 如果蛋白是正电，配体应该是负电
   - 如果蛋白是负电，配体应该是正电
3. **设计修改**
   - 添加/移除带电基团
   - 改变取代基电负性
   - 调整氢键供体/受体

---

## 输出文件说明

EC 分析会生成以下文件：

| 文件 | 说明 |
|------|------|
| `ec_map.pdb` | EC 值映射到 B-factor 的伪 PDB 文件，用于 PyMOL 可视化 |
| `ec_data.csv` | 包含所有表面点的坐标、EC 值、蛋白/配体势能 |
| `protein.pqr` | PDB2PQR 生成的蛋白 PQR 文件 |
| `protein_pot.dx` | APBS 生成的蛋白静电势网格 |

### PyMOL 可视化

EC map 在 PyMOL 中的颜色编码：
- **绿色**：正 EC（静电互补）
- **红色**：负 EC（静电冲突）
- **白色**：EC ≈ 0

#### 方法 1：离散点云可视化（默认）

```python
# 在 PyMOL 中手动加载和着色
cmd.load('ec_map.pdb', 'ec_map')
cmd.spectrum('b', 'red_white_green', 'ec_map', minimum=-100, maximum=100)
cmd.show('spheres', 'ec_map')
cmd.set('sphere_scale', 0.15, 'ec_map')
```

#### 方法 2：光滑分子表面可视化（推荐，类似文献图片）

这种方法生成类似文献中常见的光滑彩色分子表面，效果如下图所示：

![EC Surface Example](../assets/ec_surface_example.png)

```python
from gluetk.ligand_ec_calculator import calculate_ligand_ec, visualize_ec_surface

# 1. 先运行 EC 分析
result = calculate_ligand_ec('complex', 'LIG', output_dir='./ec_output')

# 2. 创建光滑表面可视化
visualize_ec_surface(
    'complex',           # PyMOL 对象名
    'LIG',               # 配体残基名
    ec_result=result,    # EC 分析结果
    surface_type='gaussian',  # 表面类型：'gaussian', 'solvent', 'molecular'
    transparency=0.0     # 透明度：0.0=不透明, 1.0=完全透明
)

# 3. 高质量渲染
cmd.ray(1920, 1080)
cmd.png('ec_surface.png', dpi=300)
```

**表面类型说明：**

| 类型 | 说明 | 适用场景 |
|------|------|----------|
| `gaussian` | 高斯表面，最光滑 | 发表级图片 |
| `solvent` | 溶剂可及表面 | 展示溶剂接触区域 |
| `molecular` | 分子表面（vdW） | 展示分子形状 |

#### 方法 3：CGO 点云表面

创建密集的彩色点云，近似光滑表面：

```python
from gluetk.ligand_ec_calculator import visualize_ec_surface_cgo

visualize_ec_surface_cgo(
    'complex',
    'LIG',
    ec_values=result['ec_values'],
    surface_points=result['surface_points'],
    point_size=0.3  # 点大小
)
```

#### 可视化效果对比

| 方法 | 效果 | 性能 | 适用场景 |
|------|------|------|----------|
| 离散球体 | ⭐⭐ | 快 | 快速预览 |
| 光滑表面 | ⭐⭐⭐⭐ | 中 | 发表图片 |
| CGO 点云 | ⭐⭐⭐ | 慢 | 大量表面点 |

#### 高质量渲染设置

```python
# 设置渲染参数
cmd.set('ray_shadow', 0)        # 关闭阴影（更清晰）
cmd.set('antialias', 2)         # 抗锯齿
cmd.set('ambient', 0.4)         # 环境光
cmd.set('spec_reflect', 0.5)    # 镜面反射
cmd.set('surface_quality', 1)   # 表面质量

# 渲染
cmd.ray(1920, 1080)
cmd.png('ec_surface_hq.png', dpi=300)
```

#### 颜色方案自定义

```python
# 使用不同的颜色方案
# 红-白-蓝（经典静电势配色）
cmd.spectrum('b', 'red_white_blue', 'ec_surface_LIG', minimum=-100, maximum=100)

# 蓝-白-红（反转）
cmd.spectrum('b', 'blue_white_red', 'ec_surface_LIG', minimum=-100, maximum=100)

# 自定义颜色
cmd.set_color('ec_negative', [1.0, 0.2, 0.2])  # 红色
cmd.set_color('ec_neutral', [1.0, 1.0, 1.0])   # 白色
cmd.set_color('ec_positive', [0.2, 0.8, 0.2])  # 绿色
```

---

## 参考文献

1. Clark, T. et al. (2007). σ-Holes. J. Mol. Model. 13, 291-296.
2. Baker, N.A. et al. (2001). Electrostatics of nanosystems. PNAS 98, 10037-10041.
3. Gasteiger, J. & Marsili, M. (1980). Iterative partial equalization of orbital electronegativity. Tetrahedron 36, 3219-3228.
4. Bauer, M.R. & Mackey, M.D. (2019). Electrostatic Complementarity as a Fast and Effective Tool to Optimize Binding and Selectivity of Protein–Ligand Complexes. J. Med. Chem. 62, 3036-3050.