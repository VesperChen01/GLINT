# 分子胶设计工作流
## Molecular Glue Design & PROTAC Development Workflow

> **理性药物设计指南** - 从 Ternary Complex 建模到分子胶优化  
> **Rational Drug Design Guide** - From Ternary Complex Modeling to Glue Optimization

---

## 🎯 工作流概览

```
新底物蛋白 (Target Protein with G-motif)
    ↓
Step 1: G-loop 对齐建模 (Alignment & Modeling)
    ↓
Step 2: 碰撞与界面分析 (Clash & Interface Analysis)
    ↓
Step 3: Exit Vector 识别 (Exit Vector Identification)
    ↓
Step 4: 静电环境分析 (Electrostatic Environment)
    ↓
Step 5: 分子胶优化建议 (Glue Optimization Strategy)
    ↓
输出: 设计报告 + PROTAC Linker 连接点
```

---

## 📦 第一步：Ternary Complex 建模

### 目标
将新底物的 G-motif 对齐到已知 CRBN-IMiD-neosubstrate 模板结构

### 常用模板

| 蛋白 | PDB | G-loop链 | G-loop范围 | CRBN链 | 配体 | 用途 |
|------|-----|----------|------------|--------|------|------|
| GSPT1 | 6H0G | A | 60-67 | E | CC9 | 通用模板 |
| CK1α | 5FQD | A | 36-43 | B | LEN | Lenalidomide |
| IKZF3 | 6H0F | A | 143-150 | E | CC9 | ZF 蛋白 |

### 使用方法

```python
# 在 PyMOL 中

# 1. 加载模板和目标结构
fetch 6H0G, async_=0        # 模板 (GSPT1)
load my_target.pdb           # 你的目标蛋白

# 2. 对齐 G-loop
result = align_gloop_for_modeling(
    template_obj='6H0G',
    template_chain='A',
    template_gloop_range='60-67',
    target_obj='my_target',
    target_chain='A',
    target_gloop_range='100-107',
    output_obj='my_target_aligned'
)

print(f"对齐 RMSD: {result['rmsd']:.3f} Å")
```

### 输出解读

- **RMSD < 1.0 Å**: ✅ 优秀，G-loop 高度相似
- **RMSD 1.0-2.0 Å**: ⚠️ 一般，可能需要侧链调整
- **RMSD > 2.0 Å**: ❌ 较差，考虑使用其他模板或手动优化

---

## 🔍 第二步：碰撞与界面分析

### 目标
检测新底物 G-loop 与 CRBN 是否有严重碰撞（clashes）

### 使用方法

```python
# 检测碰撞
clash_result = detect_clashes_at_interface(
    aligned_obj='my_target_aligned',
    aligned_chain='A',
    aligned_gloop_range='100-107',
    crbn_obj='6H0G',              # 模板中的 CRBN
    crbn_chain='E',
    clash_threshold=2.0,           # 碰撞阈值
    contact_threshold=4.5          # 接触阈值
)

print(f"碰撞数: {clash_result['clash_count']}")
print(f"严重碰撞: {clash_result['has_severe_clashes']}")
```

### 结果判断

#### 碰撞分类

| 碰撞数 | 严重碰撞 (< 1.5Å) | 判断 | 建议 |
|--------|-------------------|------|------|
| 0-2 | ❌ | ✅ 优秀 | 继续下一步 |
| 3-5 | ❌ | ⚠️ 可接受 | 考虑侧链优化 |
| > 5 | ❌ | ⚠️ 多碰撞 | 需要优化或换模板 |
| 任意 | ✅ | ❌ 不可行 | 必须解决 |

#### 解决策略

1. **轻微碰撞 (2.0-1.5 Å)**:
   - 侧链旋转异构优化
   - PyMOL: `sculpt my_target_aligned and resi 100-107`

2. **中度碰撞 (1.5-1.0 Å)**:
   - Molecular dynamics (MD) relaxation
   - 或手动调整侧链

3. **严重碰撞 (< 1.0 Å)**:
   - 可能需要不同模板
   - 或该 G-motif 不适合当前分子胶

---

## 🚀 第三步：Exit Vector 识别

### 目标
找到分子胶上适合连接 PROTAC linker 的位置

### 原理

**Exit vector** 是配体上满足以下条件的原子：
1. 远离 CRBN-G-loop 核心界面
2. 溶剂暴露（可及表面积大）
3. 化学修饰不破坏结合

### 使用方法

```python
# 识别 exit vectors
exit_vectors = identify_exit_vectors(
    ligand_obj='6H0G',
    ligand_resname='CC9',          # 配体名称（IMiD 衍生物）
    crbn_obj='6H0G',
    crbn_chain='E',
    gloop_obj='my_target_aligned',
    gloop_chain='A',
    gloop_range='100-107',
    solvent_radius=5.0
)

# 查看推荐
for ev in exit_vectors[:3]:
    print(f"原子: {ev['atom_name']}, 评分: {ev['score']}, 原因: {ev['reasons']}")
```

### 评分系统

| 评分 | 等级 | 说明 |
|------|------|------|
| > 8.0 | ⭐⭐⭐ | 优秀 exit vector |
| 6.0-8.0 | ⭐⭐ | 良好 |
| 4.0-6.0 | ⭐ | 可用 |
| < 4.0 | ❌ | 不推荐 |

### 可视化

```python
# 高亮 exit vector
for ev in exit_vectors[:1]:
    cmd.select(f"exit_vector", f"6H0G and resn CC9 and name {ev['atom_name']}")
    cmd.show('spheres', 'exit_vector')
    cmd.color('magenta', 'exit_vector')
```

---

## ⚡ 第四步：静电环境分析

### 目标
了解 G-motif 周围的电荷分布，指导分子胶化学修饰

### 使用方法

```python
# 分析静电环境
electro = analyze_electrostatic_environment(
    obj_name='my_target_aligned',
    gloop_chain='A',
    gloop_range='100-107',
    crbn_chain='E',
    sampling_radius=10.0    # 采样半径（Å）
)

print(f"净电荷: {electro['net_charge']:+d}")
print(f"疏水比例: {electro['hydrophobic_ratio']:.1%}")
print(f"建议: {electro['recommendation']}")
```

### 设计建议对照表

| 净电荷 | 疏水比例 | 环境 | 推荐修饰 | 避免 |
|--------|----------|------|----------|------|
| > +2 | <50% | 正电 | 羧酸、磺酸、磷酸 | 胺、胍 |
| < -2 | <50% | 负电 | 胺、胍、吡啶 | 羧酸 |
| -2 ~ +2 | >50% | 疏水 | 芳香、烷基、卤素 | 极性基团 |
| -2 ~ +2 | <50% | 中性 | 平衡疏水/亲水 | 过度极性 |

### 化学基团推荐

#### 正电环境（净电荷 > +2）
```
负电基团:
-COOH (羧酸)       ← 最常用
-SO₃H (磺酸)       ← 强酸性
-PO₃H₂ (磷酸)      ← 中等
-CF₃ (三氟甲基)    ← 弱吸电子
```

#### 负电环境（净电荷 < -2）
```
正电基团:
-NH₂ (伯胺)        ← 基本
-NHR (仲胺)        ← 可调
-NR₂ (叔胺)        ← 弱碱性
-C(=NH)NH₂ (胍)   ← 强碱性
```

#### 疏水环境（疏水比 > 50%）
```
疏水基团:
-Ph (苯基)         ← 芳香
-Cy (环己基)       ← 脂肪环
-CF₃ (三氟)        ← 卤代
-t-Bu (叔丁基)     ← 位阻
```

---

## 🧪 第五步：综合分析（一键工作流）

### 使用方法

```python
# 完整工作流
result = comprehensive_glue_design_analysis(
    template_obj='6H0G',
    template_crbn_chain='E',
    template_gloop_chain='A',
    template_gloop_range='60-67',
    target_obj='my_target',
    target_chain='A',
    target_gloop_range='100-107',
    ligand_resname='CC9',
    output_report='my_design_report.txt'
)

# 自动生成报告文件
```

### 报告内容

```
================================================================================
分子胶设计分析报告
Molecular Glue Design Analysis Report
================================================================================

模板: 6H0G (A:60-67)
目标: my_target (A:100-107)
配体: CC9

--------------------------------------------------------------------------------
1. G-loop 对齐
--------------------------------------------------------------------------------
RMSD: 0.850 Å
对齐质量: ✅ 优秀 (<1Å)

--------------------------------------------------------------------------------
2. 碰撞分析
--------------------------------------------------------------------------------
碰撞数: 2
接触数: 35
严重碰撞: ✅ 否

主要碰撞:
  LEU 102 ↔ ASN 351: 1.85 Å
  VAL 105 ↔ HIS 357: 1.92 Å

--------------------------------------------------------------------------------
3. G-motif 验证
--------------------------------------------------------------------------------
中心 Gly: ✅
α-turn 氢键: ✅ (2.8 Å)
RMSD: 0.85 Å
有效几何: ✅

--------------------------------------------------------------------------------
4. Exit Vector 推荐
--------------------------------------------------------------------------------
1. C14 (C) - 评分 8.5
   原因: 远离CRBN, 远离G-loop, 溶剂暴露, 碳原子
2. C15 (C) - 评分 7.0
   原因: 远离G-loop, 溶剂暴露, 碳原子
3. C10 (C) - 评分 6.5
   原因: 溶剂暴露, 碳原子

--------------------------------------------------------------------------------
5. 静电环境与设计建议
--------------------------------------------------------------------------------
净电荷: -3
疏水比例: 45%

💡 建议: 环境偏负电：考虑引入正电基团（胺、胍）

================================================================================
分析完成
================================================================================
```

---

## 💊 第六步：分子胶优化策略

### A. 基于碰撞结果的优化

#### 情况 1: 有碰撞但无严重碰撞
**策略**: 体积调整
```
碰撞残基: LEU 102 (目标) ↔ ASN 351 (CRBN)

分子胶修饰方向:
1. 在 exit vector 位置引入更小的取代基
2. 考虑使用更扁平的芳环（如噻唑代替苯）
3. 在碰撞侧引入柔性 linker
```

#### 情况 2: 严重碰撞
**策略**: 换用不同骨架或模板
```
⚠️ 当前 IMiD 骨架不适合此 G-motif

备选方案:
1. 尝试其他 E3 配体（VHL、IAP）
2. 使用不同的 CRBN 配体（如 Pomalidomide vs Lenalidomide）
3. 考虑该蛋白不适合做分子胶降解
```

### B. 基于 Exit Vector 的 PROTAC 设计

#### PROTAC Linker 连接策略

```python
# 假设 Exit Vector 是 C14
# 推荐 Linker 类型:

1. PEG Linker (亲水)
   C14 -- (OCH₂CH₂)ₙ -- E3 Ligand
   n = 2-6

2. Alkyl Linker (疏水)
   C14 -- (CH₂)ₙ -- E3 Ligand
   n = 4-10

3. Piperazine Linker (刚性)
   C14 -- Piperazine -- (CH₂)ₙ -- E3 Ligand
```

#### Linker 长度选择

| CRBN-G-loop 距离 | Linker 长度 | 推荐类型 |
|------------------|-------------|----------|
| 15-20 Å | 短 (n=2-4) | Alkyl |
| 20-30 Å | 中 (n=4-8) | PEG |
| > 30 Å | 长 (n>8) | PEG + Piperazine |

### C. 基于静电环境的化学优化

#### 实例：环境偏负电（净电荷 -3）

**当前 IMiD**: CC9 (Lenalidomide 衍生物)

**优化方向**:
1. 在 exit vector (C14) 引入胺基：`-NH-CH₂-Ph`
2. 或引入哌啶：`-piperidine`
3. 避免羧酸或磺酸

**预期效果**:
- 增强与负电残基（Asp/Glu）的盐桥
- 提高结合亲和力
- 可能提高选择性（相对 GSPT1）

---

## 📊 实战案例：优化新底物的分子胶

### 背景
- 目标蛋白：MyProtein（含 G-motif at 100-107）
- 现有 IMiD: Lenalidomide（对 MyProtein 弱效 DC50 > 10 μM）
- 目标：设计新分子胶使 DC50 < 1 μM

### Step-by-Step 流程

```python
# 1. 建模
fetch 6H0G, async_=0
load MyProtein.pdb

align_result = align_gloop_for_modeling(
    '6H0G', 'A', '60-67',
    'MyProtein', 'A', '100-107'
)
# RMSD: 1.2 Å (可接受)

# 2. 碰撞检测
clash_result = detect_clashes_at_interface(
    'MyProtein_aligned', 'A', '100-107',
    '6H0G', 'E'
)
# 碰撞: 4 处，严重碰撞: 无
# 主要碰撞: PHE 103 ↔ HIS 357 (1.9 Å)

# 3. Exit vector
exit_vectors = identify_exit_vectors(
    '6H0G', 'LEN', '6H0G', 'E',
    'MyProtein_aligned', 'A', '100-107'
)
# 推荐: C5 (评分 8.0) - 远离界面，溶剂暴露

# 4. 静电环境
electro = analyze_electrostatic_environment(
    'MyProtein_aligned', 'A', '100-107', 'E'
)
# 净电荷: +4, 疏水比: 35%
# 建议: 环境偏正电，引入负电基团
```

### 设计方案

#### 方案 A: 减小体积（解决碰撞）
```
原 Lenalidomide: 苯环在 exit vector 侧
↓
新设计: 替换为噻唑（更小）
预期: 减少 PHE 103 碰撞，提高亲和力 2-5 倍
```

#### 方案 B: 静电互补
```
在 C5 引入羧酸: -COOH
↓
与 MyProtein G-loop 周围的 Arg/Lys 形成盐桥
预期: 提高选择性，DC50 < 5 μM
```

#### 方案 C: PROTAC 化（最终方案）
```
C5-羧酸修饰 → PEG4 linker → VHL ligand
↓
Heterobifunctional PROTAC
预期: DC50 < 100 nM
```

---

## 🔬 实验验证建议

### 体外实验（In Vitro）

1. **SPR/ITC**: 测量分子胶与 CRBN 亲和力
2. **AlphaLISA**: 检测 ternary complex 形成
3. **Western Blot**: 验证底物降解（DC50、Dmax）

### 细胞实验（Cell-Based）

1. **Degradation Assay**: 
   - 时间依赖性（0-24h）
   - 浓度依赖性（0.1 nM - 10 μM）
   
2. **Selectivity Panel**:
   - 对照已知底物（GSPT1, CK1α）
   - 确认选择性窗口

3. **Rescue Experiment**:
   - + MLN4924 (抑制 Cullin neddylation)
   - + CRBN CRISPR KO
   - 确认 Cullin-RING 依赖性

---

## 📖 参考文献 & 工具

### 关键文献
1. Ternary complex 结构优化：Gadd, M.S. et al. (2017) *Nat. Chem. Biol.* 13, 514–521
2. Exit vector 设计：Nowak, R.P. et al. (2018) *Nat. Chem. Biol.* 14, 706–714
3. PROTAC linker 优化：Crew, A.P. et al. (2018) *J. Med. Chem.* 61, 583–598

### 在线工具
- **PROTAC-DB**: http://cadd.zju.edu.cn/protacdb/
- **Degrader Linker Library**: Tocris Bioscience
- **CRBN Structure**: PDB 6H0G, 5FQD, 6H0F

---

## 💡 最佳实践

1. **先验证，后设计**: 确认 G-motif 是标准结构再优化分子胶
2. **多模板对照**: 尝试不同模板（GSPT1/CK1α/IKZF3）找最佳匹配
3. **迭代优化**: 碰撞 → Exit vector → 静电 → 合成 → 测试 → 再优化
4. **平行设计**: 同时设计 2-3 个方案，降低风险
5. **文献先行**: 查看该蛋白是否已有分子胶/PROTAC 报道

---

**版本**: 1.0  
**作者**: Vesper (GlueTK)  
**最后更新**: 2025-11-11  
**适用**: 分子胶降解剂 & PROTAC 设计
