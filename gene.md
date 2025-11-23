# 项目标题

基于公开 API 的「疾病 → 关联蛋白 → 分子胶潜力」分析，并在 PyMOL 中可视化(新增一个页面模块)
，方便结构研究者按疾病筛选潜在分子胶靶点。

---

## 一、整体目标（给工程师看的高层描述）

做一个 **Python + PyMOL 插件**，满足下面这个使用体验：

1. 用户在 PyMOL 里输入一条命令（或在插件界面里填）：

   * 疾病名称（英文，例如：`"multiple myeloma"` / `"acute myeloid leukemia"`）
   * 以及候选靶点数量（例如 Top 50）。
2. 插件调用公开 API（推荐：**Open Targets Platform**），自动完成：

   * 根据疾病名称搜索疾病 ID（EFO 等）；
   * 用疾病 ID 查询该疾病的 **关联靶点列表**（disease–target association）；
   * 取出前 N 个靶点及其疾病关联分数。
3. 对选中的某一个靶点，用户可以：

   * 在 PyMOL 中加载该靶点的结构（PDB 或 AlphaFold2）；
   * 调用命令，在结构上显示：

     * **疾病关联程度**（来自 Open Targets 的 association score 或我们后续附加的 DRI 分数）；
     * **与 E3（例如 CRBN）相关的兼容度指标**（如果后续加入表达/细胞系信息的话）；
     * **分子胶相关的结构热点**（残基级评分，由现有 gloop/结构分析流程生成）。
4. 插件在 PyMOL 里：

   * 在控制台打印这个靶点的综合评分；
   * 在结构上高亮残基级“胶点热点”（例如 SurfHot_score ≥ 阈值的残基，用红色球体显示）。

**总结：**

> * 接收任意疾病名；
> * 用 API 找到关联蛋白；
> * 和本地已有的“分子胶结构打分文件”结合，在结构里可视化。

---

## 二、外部数据与 API 选择

### 1. 疾病 → 靶点：Open Targets Platform API（推荐）

* 官网：Open Targets Platform
* 提供 GraphQL API，可以

  * 搜索疾病：`search(queryString: ..., entityTypes: [DISEASE])`
  * 查询疾病的 associatedTargets：`disease(efoId: ...) { associatedTargets {...} }`
* 关键信息：

  * 疾病 ID（通常是 EFO ID，例如多发性骨髓瘤 `EFO_0001378`）
  * 靶点基因 symbol（`approvedSymbol`）
  * 疾病–靶点总体关联分数 `overallAssociationScore`

**期望行为：**

* 函数 `search_disease(name: str) -> List[CandidateDisease]`

  * 输入：任意疾病名称字符串（英文）
  * 输出：若干候选疾病（ID + name），供用户选择一个。

* 函数 `get_disease_targets(disease_id: str, topN: int) -> DataFrame`

  * 输入：疾病 ID（如 EFO_0001378）、topN
  * 输出：至少包含以下列的 DataFrame：

    * `symbol`（基因符号）
    * `overallAssociationScore`（疾病–靶点关联分数）

> 说明：后续如果需要可加入更多证据维度（expression, genetic, drug 等），但 V1 只需 association 总分即可。

### 2. 分子胶 / 结构侧数据（已存在，由现有流程提供）

用户已有独立的结构分析流程（如 gloop），会为每个蛋白生成残基级的“胶点评分”。我们假设：

* 对每个靶点基因（symbol 或 UniProt ID），已有一个本地 CSV 文件：

  * 文件命名示例：`residue_scores_{SYMBOL}.csv` 或 `residue_scores_{UNIPROT}.csv`
  * 主要列：

    * `residue_index`：残基序号（与 PDB/AF2 中的 resi 对应）
    * `SurfHot_score`：该残基作为分子胶热点的评分（0~1）
    * （可选）`motif_flag`、`glue_patch_score` 等

* 对应的结构文件（PDB/AF2）：

  * 命名示例：`AF_{UNIPROT}.pdb`
  * 用户会在 PyMOL 中自行 `load`，并给一个对象名（例如与 symbol 一致）。

**工程需求：**

* 插件要有办法：

  * 根据基因 symbol 找到对应的 `residue_scores_*.csv` 文件；
  * 读出 `SurfHot_score`，并根据阈值在 PyMOL 对象上高亮。

### 3. E3 / 细胞系相关信息（可选，预留接口）

短期内可以只用 Open Targets 的 `overallAssociationScore` 作为疾病相关性指标。

未来版本预留接口：

* 用另一个脚本/模块，根据基因 symbol 计算：

  * 与 E3（例如 CRBN）的共表达/共存在得分 `E3_score`；
* 输出到一个本地 CSV，例如：`gene_level_scores.csv`：

  * 列：`symbol`, `E3_score`, `extra_disease_score` 等

插件需要能够：

* 若该表存在，则在 PyMOL 输出中显示 `E3_score`；
* 若不存在，则跳过该项，仅显示 Open Targets 的分数。

---

## 三、功能模块设计（工程师视角）

### 模块 1：疾病搜索与选择

**功能：** 将疾病名称字符串转换为具体疾病 ID。

* 输入：疾病名称（英文字符串）

* 步骤：

  1. 调用 Open Targets 的 `search` 接口：

     * `entityTypes: [DISEASE]`
     * `queryString: 输入疾病名`
  2. 返回候选疾病列表（ID + name），按相关性排序。
  3. 若结果不唯一：

     * 在命令行打印候选列表，并提示用户选择一个 ID；
     * 或在 GUI 中拉一个下拉框供选择。

* 输出：选定的疾病 ID（例如 `EFO_0001378`）

### 模块 2：疾病 → 靶点列表

**功能：** 从 Open Targets 获取指定疾病的关联靶点列表。

* 输入：疾病 ID（EFO）、Top N

* 步骤：

  1. 调用 GraphQL：`disease(efoId: ...) { associatedTargets(...) }`
  2. 解析返回结果，提取：

     * `target.approvedSymbol` 作为基因 symbol
     * `overallAssociationScore` 作为疾病–靶点关联分数
  3. 生成 pandas DataFrame / Python 列表。

* 输出：

  * `targets_df`，列：`symbol`, `overallAssociationScore`

* 额外：将结果保存在本地 CSV，例如：

  * `targets_{disease_id}.csv`

### 模块 3：PyMOL 插件命令（疾病层）

在 PyMOL 中注册一个命令，例如：

```text
ot_disease_targets "multiple myeloma", 30
```

含义：

* 使用 Open Targets 搜索 `"multiple myeloma"` 对应的疾病 ID；
* 提取该疾病 Top 30 个关联靶点；
* 在 PyMOL 控制台打印一个表格，如：

```text
Rank  Symbol   AssocScore
1     IKZF1    0.82
2     GSPT1    0.78
...
```

* 同时把这个列表保存到 CSV 文件，方便用户离线查看。

### 模块 4：PyMOL 插件命令（单靶点 + 结构可视化）

在 PyMOL 注册第二个命令，例如：

```text
ot_glue_insight protein_obj, gene_symbol, threshold=0.7
```

功能：

1. 在本地数据中查找该基因的疾病相关性信息：

   * 第一优先：刚刚从 Open Targets 获取的 `targets_{disease_id}.csv` 中是否存在该基因；
   * 若存在，读取其 `overallAssociationScore`；
   * 若预留扩展表 `gene_level_scores.csv` 存在，则一并读取并显示 `E3_score` 等。
2. 在本地路径查找残基级打分文件：

   * 按约定命名：`residue_scores_{gene_symbol}.csv`
   * 若找到：

     * 读取 `residue_index`, `SurfHot_score`
     * 在给定的 PyMOL 对象 `protein_obj` 上：

       * 将所有残基先统一设置为基础颜色（如浅灰）
       * 对 `SurfHot_score ≥ threshold` 的残基：

         * 建立 selection，例如：`{protein_obj}_{gene_symbol}_hot`
         * 染色为红色 `color red, selection`
         * 显示为球体 `show spheres, selection`
3. 在 PyMOL 控制台打印汇总信息，例如：

```text
[ot_glue_insight]
Gene: IKZF1
Disease assoc score (OpenTargets): 0.82
E3_score (if available): 0.75
Hotspot threshold: 0.70
#hot residues: 15
```

---

## 四、目录与配置建议

建议工程师按如下方式组织项目（仅示意）：

```text
project_root/
  config.py                  # API endpoint、token、数据路径等配置
  open_targets_api.py        # 与 Open Targets 交互的封装
  gene_level_scores.csv      # （可选）疾病+E3 等本地扩展分数
  residue_scores/            # 残基级打分CSV目录
    residue_scores_IKZF1.csv
    residue_scores_GSPT1.csv
    ...
  pymol_plugin/
    __init__.py
    ot_disease_targets.py    # PyMOL命令：疾病→靶点列表
    ot_glue_insight.py       # PyMOL命令：单靶点结构可视化
```

`config.py` 示例内容：

```python
OPEN_TARGETS_GRAPHQL_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
RESIDUE_SCORES_DIR = "/absolute/path/to/residue_scores"
GENE_LEVEL_SCORES_CSV = "/absolute/path/to/gene_level_scores.csv"  # 可选
```

---

## 五、PyMOL 集成方式（工程说明）

* 插件本质上是若干 Python 函数 + `cmd.extend` 调用。
* 示例：

```python
from pymol import cmd

from open_targets_api import search_disease, get_disease_targets
from config import RESIDUE_SCORES_DIR

# 疾病→靶点列表

def ot_disease_targets(disease_name, topN=30):
    # 1. 调用 search_disease，得到疾病ID
    # 2. 调用 get_disease_targets，得到DataFrame
    # 3. 打印TopN表格；保存CSV
    ...

cmd.extend("ot_disease_targets", ot_disease_targets)


# 单靶点 + 结构可视化

def ot_glue_insight(protein_obj, gene_symbol, threshold=0.7):
    # 1. 从 gene_level_scores / OpenTargets CSV 找该gene的疾病相关信息
    # 2. 从 RESIDUE_SCORES_DIR 中加载 residue_scores_{gene_symbol}.csv
    # 3. 在 protein_obj 上根据 SurfHot_score ≥ threshold 染色
    ...

cmd.extend("ot_glue_insight", ot_glue_insight)
```

使用示例（用户角度）：

```pymol
run /path/to/pymol_plugin/ot_disease_targets.py
run /path/to/pymol_plugin/ot_glue_insight.py

# 1) 搜索疾病的关联靶点
ot_disease_targets "multiple myeloma", 30

# 2) 用户根据输出列表，选择一个靶点（例如 IKZF1），手动加载结构
load AF_Q13422.pdb, IKZF1

# 3) 对该靶点做分子胶热点可视化
ot_glue_insight IKZF1, IKZF1, 0.8
```

---

## 六、版本划分建议（迭代计划）

### V1（最小可用版本）

* 支持：

  * 用 Open Targets 搜索疾病并获取 TopN 靶点（symbol + associationScore）；
  * PyMOL 中通过 `ot_disease_targets` 命令打印结果；
  * 读取本地 `residue_scores_{symbol}.csv`，在结构中按阈值染色热点；
  * 控制台输出：疾病关联分数（来自 Open Targets 的 overallAssociationScore）。
* 不做：

  * 复杂的 E3_score 计算；
  * 多数据源整合；
  * GUI 面板。

### V2（增强版）

* 增加：

  * 读取 `gene_level_scores.csv`，显示额外的疾病相关性/CRBN 兼容度分数（E3_score）。
  * 按综合分数（例如 0.6*OpenTargets + 0.4*E3_score）对靶点再排序。

### V3（可选高级版）

* 增加一个简单 GUI：

  * 在 PyMOL 右侧面板中选择疾病、修改 TopN、调整阈值；
  * 点击按钮进行查询和可视化，而不是只靠命令行。

---

## 七、给工程师的备注

1. 用户不要求 Web 界面，只要本地 PyMOL 环境能跑即可。
2. 重点是 **“疾病可搜索，不写死多发性骨髓瘤”**，所有逻辑围绕“疾病名→疾病ID→靶点列表”展开。
3. 分子胶结构打分（残基级 SurfHot_score 等）由用户现有流程生成，工程侧只负责：

   * 按约定命名和路径读取这些 CSV；
   * 在 PDB/AF2 结构上正确对应 residue index 并上色。
4. 若需要示例数据（假 CSV）来开发/调试，用户可以提供一个小的 `residue_scores_IKZF1.csv` 和一个对应的 AF2/PDB 文件，方便先实现插件逻辑。