# GlueTK：分子胶降解剂分析的计算框架（中文对照版 Outline）

**目标期刊**：Journal of Chemical Information and Modeling (JCIM)  
**文章类型**：Application Note / Software  
**预计长度**：≤5000 words（JCIM Application Note 规则，含图表折算）  
**最后更新**：2025-11-11

---

## JCIM Application Note 合规要点（写作前先锁死）

- **标题必须包含应用名**：标题里清楚写出 `GlueTK`
- **篇幅上限**：≤5000 words（含 abstract + 正文 + 图形折算；参考文献不计入）
- **至少 1 张图**：建议保留 3–6 张关键图，但要做字数折算与取舍
- **强调“科学/技术/可用性进步”**：用 3–5 条 bullet 明确你比现有方案多了什么（整合性、可复现输出、面向设计的口袋/Neo-表位等）
- **操作系统无关（Windows/Linux/macOS）**：Methods 或 Implementation 明确写 cross-platform，并在 SI/README 给安装说明
- **可获取性**：论文发表时软件需可公开获取（开源仓库 + release/installer + 版本号 +（建议）Zenodo DOI）
- **评审可测试性**：需要能提供 reviewer 测试（按需提供下载与测试数据/脚本，同时保护 reviewer 匿名）

### 字数核算规则（写作时必须同步记账）

- **单栏图/示意图/表**：按 **300 words/张** 估算
- **双栏图/示意图/表**：按 **600 words/张** 估算
- Word 里用 “Tables” 功能做的表：**会计入 word count**；建议把复杂表改成图片（按 300/600 折算）更可控

### Word count statement 模板（投稿时要交的那段话）

> Word count statement: The total word count of this Application Note is **XXXX** words, including the abstract, main text, and graphics. Graphics were approximated using **300 words per single-column figure/table** and **600 words per double-column figure/table**, following the journal guidelines. References were excluded from the count.

---

## 题目备选（Title Options）

1. **GlueTK：基于 PyMOL 的分子胶机制分析插件——蛋白-蛋白界面、Neo-表位识别与口袋引导设计**
**推荐**：选项 1（突出 pocket analysis 新功能与“机制分析+设计启发”定位）

---

## 摘要（Abstract，≤250 words）

### 建议结构
1. **背景（2–3 句）**
   - 分子胶降解剂是 TPD（靶向蛋白降解）中快速发展的药物类型
   - 当前缺乏面向分子胶“界面诱导/Neo-表位/口袋”特征的一体化计算与可视化工具
   -（如仍保留对比叙事）现有工作流难以系统性评估三元复合物的结构依据与设计要点

2. **方法（3–4 句）**
   - 我们开发 GlueTK：一个在 PyMOL 内运行的分子胶分析插件
   - 核心功能包括：PPI 界面与 BSA、Neo-表位、CRBN G-motif/G-loop、界面口袋检测与 glue-pocket 关联、（可选）协同评分
   - 交互判据对齐发表标准（例如 Maestro 级别氢键距离阈值）

3. **结果（2–3 句）**
   - 在多个已知分子胶三元复合物上展示：界面解析、Neo-表位定位、口袋与结合位点的关联、以及可直接用于论文作图的输出
   - 报告运行时间与输出形式（CSV/图形/对象），强调可复现与易用性

4. **结论（1–2 句）**
   - GlueTK 提供一个可复现、可视化、面向设计启发的分子胶结构分析工作流
   - 开源发布，附完整文档与 GUI

---

## 1. 引言（Introduction，800–1000 words）

### 1.1 TPD 与分子胶概述
- TPD 的临床与产业背景
- 分子胶与 PROTAC 的差异（可简要带过；若主线改为“工作流”，不必强调二分类）
- 代表性分子胶体系（CRBN、DCAF15 等）

### 1.2 结构机制要点（分子胶视角）
- 分子胶的关键结构特征：诱导/稳定 PPI、形成/增强界面、产生 Neo-表位、界面口袋与形状互补
- 设计关注点：哪些残基/口袋可能决定选择性与 SAR

### 1.3 现有计算工具与缺口
- 现有工具多偏向通用蛋白-配体相互作用或数据库检索
- 分子胶特异的整合缺口：
  - PPI 强度与界面相互作用的统一输出
  - Neo-表位（底物侧）识别与可视化
  -（CRBN）G-motif/G-loop 这类 E3 特异招募特征
  - 界面口袋、glue-pocket 关联与可用于优化的指标

### 1.4 本工作的贡献（Our Contribution）
- GlueTK：整合 PPI + Neo-表位 +（可选）G-motif + pocket analysis 的一站式 PyMOL 工作流
- 输出面向论文与 SI：标准化阈值、CSV 报告、可直接截图的 PyMOL 对象与配色
- 跨平台、开源、GUI+命令行双入口

---

## 2. 方法（Methods，1500–2000 words）

### 2.1 算法设计（Algorithm Design）

#### 2.1.1 蛋白-蛋白界面（PPI）检测
**算法 1：PPI 界面分析**
- 输入：结构（E3 与底物链），距离阈值（默认 4.5 Å）
- 输出：界面残基对、相互作用类型统计、BSA（可计算时）、界面强度评分
- 关键参数：界面距离阈值、强界面阈值（例如 BSA > 800 Å² 或接触数/强度综合）

#### 2.1.2 Neo-表位识别
**算法 2：Neo-底物表位检测**
- 核心逻辑：从 glue 的“桥接原子”出发，筛选同时与 glue 与 E3 接触的底物残基
- 输出：Neo-表位残基列表、桥接原子数、置信度（0–1）
- 建议在 SI 给出阈值敏感性：距离阈值 4–6 Å 的稳定性范围

#### 2.1.3 CRBN G-motif / G-loop（可选模块）
**算法 3：G-motif 识别（结构模板 + RMSD）**
- 输入：底物链与 loop 区段、模板模式（builtin/ideal/custom）
- 输出：匹配结果、RMSD、对应残基
- 生物意义：帮助解释 CRBN 招募特异性，并为“可招募底物筛选”提供线索

#### 2.1.4 口袋检测与 glue-pocket 关联
**算法 4：界面口袋分析**
- 输入：结构、界面区域、glue 原子
- 输出：口袋集合（体积/深度/疏水性等）、与 glue 的距离/重叠指标、关键口袋残基
- 应用：给出可用于结构优化的“口袋-残基列表”和可视化

#### 2.1.5（可选）协同/评分模型
- 若保留：明确它是“启发式/经验性”评分，用于排序与摘要，不作为严格自由能
- 在 Discussion 里承认局限与未来方向（MM/PBSA、FEP、ML potential 等）

### 2.2 实现细节（Implementation Details）

#### 2.2.1 软件架构
- 语言：Python（PyMOL 插件）
- 依赖：PyMOL、NumPy/SciPy/Matplotlib、RDKit（可选）、Vina（可选）
- GUI：PyQt5/6
- 模块化：各分析模块可独立调用（命令/GUI）

#### 2.2.2 相互作用判据（可发表标准）
**表 2：相互作用参数（示例）**
- 氢键距离/角度阈值
- 盐桥、疏水、π-π 等判据

### 2.3 验证数据集（Validation Dataset）
（如果主线为“工作流”，建议标题改为“Case Studies & Reproducibility Set”）
- 选取 2–4 个代表性分子胶三元复合物（CRBN 系 + DCAF15 系各至少一个）
- 表格字段：PDB ID、链 ID、glue resname、文献来源、验证要点（Neo-表位/关键口袋/关键相互作用）

---

## 3. 结果（Results，1000–1500 words）

### 3.1 案例 1：CRBN 体系（例如 CC-885）
- 运行流程（PyMOL 命令）：
  - `ppi_analyze`
  - `neo_epitope_find`
  -（可选）`find_crbn_g_motif` / `analyze_g_motif_glue_binding`
  - `detect_pockets` / `comprehensive_glue_pocket_analysis`
- 展示点：
  - 界面残基与关键相互作用
  - Neo-表位残基定位（与文献/结构直觉一致）
  - 口袋与 glue 的空间关联 → 提出可优化的口袋残基/区域

### 3.2 案例 2：DCAF15 体系（例如 Indisulam）
- 同一流程跑一遍，强调“跨 E3 的可迁移性”

### 3.3 可复现性与参数稳定性（轻量）
- 重跑一致性：同结构重复运行，关键输出一致
- 参数敏感性：距离阈值等在合理范围内结论稳定

### 3.4 计算性能
- 每个模块耗时、总耗时、输出文件

---

## 4. 讨论（Discussion，800–1000 words）

### 4.1 优势
- 一站式：结构载入 → 分析 → 输出 → 作图
- 面向分子胶：Neo-表位、（可选）CRBN 特征、界面口袋与设计启发
- 输出可复现（CSV）且可发表（可视化）

### 4.2 局限与未来方向
- 静态结构（未覆盖动力学/多构象）
- 评分为经验性（若保留）
- 口袋检测对复杂界面仍需人工核对
- 扩展：更多 E3 特异模块、轨迹分析、更多数据集

### 4.3 与现有工具比较（可保留，但把“二分类”条目弱化）
- 重点对比：是否支持 Neo-表位、界面口袋、PyMOL 内整合、可视化与批量输出

---

## 5. 结论（Conclusions，200–300 words）
- 总结 GlueTK 的工作流价值与可复现性
- 强调开源、可用性与面向设计的输出

---

## 支持信息（Supporting Information）
- SI-1：算法伪代码与参数说明
- SI-2：案例数据集表（PDB/链/resname/参考）
- SI-3：阈值敏感性与重跑一致性
- SI-4：安装与快速教程
- SI-5：示例输出文件（CSV/PNG）

---

## 数据与代码可用性（Data and Code Availability）
- Code：GitHub release（建议配 Zenodo DOI）
- Data：案例结构来源（PDB IDs）+ 生成的 CSV 与作图源数据
