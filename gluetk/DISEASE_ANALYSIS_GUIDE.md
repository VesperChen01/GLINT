# Open Targets 疾病靶点分析 - 使用指南

## 功能概述

GlueTK 现已集成 Open Targets Platform API，支持基于疾病的分子胶靶点发现工作流：

1. **疾病搜索** → 根据疾病名称查找 EFO ID
2. **靶点查询** → 获取疾病关联的蛋白靶点列表
3. **热点可视化** → 在 PyMOL 中显示分子胶结合热点

---

## 快速开始

### 1. 安装依赖

```bash
pip install requests pandas
```

### 2. 在 PyMOL 中使用

#### 步骤 1: 查询疾病关联靶点

```pymol
# 搜索多发性骨髓瘤的关联靶点（Top 30）
ot_disease_targets disease_name="multiple myeloma", top_n=30

# 搜索急性髓系白血病的关联靶点
ot_disease_targets disease_name="acute myeloid leukemia", top_n=50
```

**输出：**
- PyMOL 控制台显示靶点列表（基因符号、关联分数）
- 保存 CSV 文件到 `./disease_data/targets_{disease_id}.csv`

#### 步骤 2: 加载靶点结构

```pymol
# 方法 1: 从 PDB 加载
fetch 6h0f, IKZF1

# 方法 2: 加载 AlphaFold 结构
load AF_Q13422.pdb, IKZF1
```

#### 步骤 3: 可视化分子胶热点

```pymol
# 前提：需要准备残基评分文件（见下文）
ot_glue_insight protein_obj="IKZF1", gene_symbol="IKZF1", threshold=0.7
```

**效果：**
- 热点残基显示为红色球体
- 控制台输出疾病关联分数和热点统计

---

## 数据准备

### 残基评分文件格式

在使用 `ot_glue_insight` 之前，需要准备残基评分文件：

**文件路径：** `./disease_data/residue_scores/residue_scores_{SYMBOL}.csv`

**文件格式：**
```csv
residue_index,SurfHot_score,motif_flag
23,0.85,1
45,0.92,1
67,0.65,0
89,0.78,1
```

**列说明：**
- `residue_index`: 残基序号（对应 PDB/AlphaFold 中的 resi）
- `SurfHot_score`: 分子胶热点评分（0-1）
- `motif_flag`: 是否为 G-loop motif（可选）

### 生成残基评分文件

使用 GlueTK 的 G-loop 分析功能生成：

```pymol
# 分析 G-loop motif
find_crbn_g_motif obj_name="IKZF1", output_csv="./disease_data/residue_scores/residue_scores_IKZF1.csv"
```

---

## 命令参考

### `ot_disease_targets`

查询指定疾病的关联靶点列表。

**语法：**
```pymol
ot_disease_targets disease_name=<疾病名称>, top_n=<数量>, output_dir=<输出目录>
```

**参数：**
- `disease_name` (必需): 疾病名称（英文），如 "multiple myeloma"
- `top_n` (可选): 返回的靶点数量，默认 30
- `output_dir` (可选): 输出目录，默认 `./disease_data`
- `auto_select` (可选): 自动选择唯一结果，默认 True

**示例：**
```pymol
ot_disease_targets disease_name="breast cancer", top_n=50
ot_disease_targets disease_name="Alzheimer disease", top_n=20, output_dir="./my_data"
```

---

### `ot_glue_insight`

可视化单个靶点的分子胶热点。

**语法：**
```pymol
ot_glue_insight protein_obj=<对象名>, gene_symbol=<基因符号>, threshold=<阈值>
```

**参数：**
- `protein_obj` (必需): PyMOL 对象名称
- `gene_symbol` (必需): 基因符号（如 IKZF1）
- `threshold` (可选): 热点阈值，默认 0.7
- `scores_dir` (可选): 残基评分文件目录
- `disease_data_dir` (可选): 疾病数据目录

**示例：**
```pymol
ot_glue_insight protein_obj="IKZF1", gene_symbol="IKZF1", threshold=0.8
ot_glue_insight protein_obj="my_protein", gene_symbol="GSPT1", threshold=0.6
```

---

## 完整工作流示例

### 示例：多发性骨髓瘤靶点分析

```pymol
# 1. 查询疾病关联靶点
ot_disease_targets disease_name="multiple myeloma", top_n=30

# 输出示例：
# Rank  Symbol   AssocScore
# 1     IKZF1    0.8234
# 2     GSPT1    0.7845
# 3     ...

# 2. 加载 Top 1 靶点结构（IKZF1）
fetch 6h0f, IKZF1

# 3. 准备残基评分文件（使用 G-loop 分析）
find_crbn_g_motif obj_name="IKZF1", output_csv="./disease_data/residue_scores/residue_scores_IKZF1.csv"

# 4. 可视化分子胶热点
ot_glue_insight protein_obj="IKZF1", gene_symbol="IKZF1", threshold=0.7

# 输出示例：
# Gene: IKZF1
# Disease Association Score: 0.8234
# Hotspot Threshold: 0.70
# Number of Hotspot Residues: 15
# Top 5 Hotspots:
#   Residue 45: 0.920
#   Residue 23: 0.850
#   ...
```

---

## 数据文件说明

### 输出文件结构

```
./disease_data/
├── targets_EFO_0001378.csv          # 疾病靶点列表
├── cache/                            # API 查询缓存
│   └── *.json
└── residue_scores/                   # 残基评分文件
    ├── residue_scores_IKZF1.csv
    ├── residue_scores_GSPT1.csv
    └── ...
```

### 靶点列表 CSV 格式

**文件：** `targets_{disease_id}.csv`

**列：**
- `symbol`: 基因符号
- `name`: 基因全名
- `score`: 疾病-靶点关联分数（0-1）
- `uniprot_id`: UniProt ID
- `ensembl_id`: Ensembl ID

---

## 常见问题

### Q1: 如何搜索中文疾病名称？

A: Open Targets API 仅支持英文查询。请使用英文疾病名称，例如：
- "multiple myeloma" (多发性骨髓瘤)
- "acute myeloid leukemia" (急性髓系白血病)
- "breast cancer" (乳腺癌)

### Q2: 找不到残基评分文件怎么办？

A: 请先使用 GlueTK 的 G-loop 分析功能生成评分文件：

```pymol
find_crbn_g_motif obj_name="your_protein", output_csv="./disease_data/residue_scores/residue_scores_SYMBOL.csv"
```

或手动创建符合格式的 CSV 文件。

### Q3: 如何查看缓存的查询结果？

A: 查询结果自动缓存到 `./disease_data/cache/` 目录，以 JSON 格式存储。删除缓存文件可强制重新查询。

### Q4: 网络请求失败怎么办？

A: 检查网络连接，确保可以访问 `https://api.platform.opentargets.org`。如果问题持续，可以使用缓存的结果或稍后重试。

---

## 技术细节

### API 端点

- **GraphQL API:** `https://api.platform.opentargets.org/api/v4/graphql`
- **文档:** https://platform-docs.opentargets.org/

### 缓存机制

- 查询结果自动缓存到本地 JSON 文件
- 缓存键基于查询内容的哈希值
- 删除缓存文件可强制重新查询

### 可视化配置

可在 `gluetk/disease_config.py` 中修改：

```python
DEFAULT_HOTSPOT_THRESHOLD = 0.7    # 默认阈值
HOTSPOT_COLOR = "red"               # 热点颜色
HOTSPOT_REPRESENTATION = "spheres"  # 显示方式
```

---

## 下一步开发

### V2 功能（计划中）

- [ ] E3 兼容度评分（CRBN 共表达分析）
- [ ] 综合评分排序（疾病关联 + E3 兼容度）
- [ ] GUI 界面集成

### V3 功能（可选）

- [ ] 交互式疾病选择界面
- [ ] 批量靶点分析
- [ ] 导出报告功能

---

## 反馈与支持

如有问题或建议，请联系 GlueTK 开发团队。

**版本：** 1.0.0  
**更新日期：** 2025-11-22
