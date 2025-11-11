# G-loop 验证改进指南
## 基于 Annual Review of Pharmacology and Toxicology 2023

### 📚 文献依据

基于文献 "CRBN-MGD Complexes Recognize a Structural Motif for Degradation" 的关键发现，我们改进了 G-loop/G-motif 的验证功能。

---

## 🔬 关键科学发现

### 1. G-loop 结构定义

G-loop 是一个 **8 个氨基酸的片段**，具有以下特征：
- **位置命名**: G-5, G-4, G-3, G-2, G-1, **G** (不变 Gly), G+1, G+2
- **结构类型**: alpha-turn 或 beta-hairpin
- **关键特征**: G-4 和 G 位之间有 backbone 氢键

### 2. CRBN 相互作用的分子基础

#### 三个关键氢键（Schrödinger 标准）

文献明确指出 G-loop 与 CRBN 形成 **3 个保守的 backbone 氢键**：

1. **G-3 羰基 O** ↔ **CRBN Asn351 侧链** (NH2)
2. **G-2 羰基 O** ↔ **CRBN His357 侧链** (ND1/NE2)
3. **G-1 羰基 O** ↔ **CRBN Trp400 侧链** (NE1)

**深度突变扫描研究**证实：
- Asn351, His357, Trp400 是获得性 MGD 耐药的热点突变位点
- CRBN N351D 突变体（失去氢键供体能力）不再降解 GSPT1

#### 侧链相互作用

文献指出在 **G-4, G-3, G-2, G+1** 位的侧链与 CRBN 形成：
- 有利的 van der Waals (vdW) 接触
- 部分情况下的极性接触

#### 保守 Glycine 的作用

- **保守 Gly 面向 MGD 药物**
- 与 IMiD (Immunomodulatory imide drugs) 形成 vdW 接触
- G 位突变为其他氨基酸会显著降低结合

---

## 🛠️ 改进功能

### 新增函数：`validate_crbn_hbonds()`

验证 G-loop 是否符合标准 CRBN 结合模式。

#### 参数
```python
validate_crbn_hbonds(
    obj_name,              # PyMOL 对象名
    g_motif_chain,         # G-loop 所在链
    g_motif_resi_range,    # 残基范围 "60-67" 或 (60, 67)
    crbn_chain,            # CRBN 链 ID
    max_hbond_dist=3.5,    # 氢键距离阈值（Å）
    min_donor_angle=120.0  # 预留角度参数
)
```

#### 返回值
```python
{
    'N351_hbond': {'found': bool, 'distance': float, 'g_pos': -3},
    'H357_hbond': {'found': bool, 'distance': float, 'g_pos': -2},
    'W400_hbond': {'found': bool, 'distance': float, 'g_pos': -1},
    'total_hbonds': int,              # 0-3
    'is_canonical_gloop': bool        # 至少 2/3 氢键
}
```

#### 示例输出
```
======================================================================
CRBN G-loop 氢键验证 (CRBN-G-loop H-bond Validation)
======================================================================
✅ G-3 → Asn351         距离: 2.75 Å
✅ G-2 → His357         距离: 2.82 Å
✅ G-1 → Trp400         距离: 2.68 Å

总氢键数: 3/3
标准 G-loop: ✅ 是
======================================================================
```

### 改进函数：`analyze_g_motif_glue_binding()`

新增参数 `validate_hbonds=True`，整合氢键验证。

#### 新增检测内容

1. **CRBN 氢键验证**（如启用）
   - 调用 `validate_crbn_hbonds()`
   - 如果 < 2 个氢键，发出警告

2. **Gly 位与 MGD 的 vdW 接触**
   - 检测保守 Gly（第 6 位）与 Glue 药物的接触
   - vdW 阈值：4.5 Å
   - 记录最小距离

3. **侧链位点与 CRBN 接触**
   - 检测 G-4, G-3, G-2, G+1 侧链与 CRBN 的接触
   - 排除 backbone 原子 (N, CA, C, O)
   - 阈值：5.0 Å（可调）

#### 返回值新增字段
```python
{
    # ... 原有字段 ...
    'crbn_hbond_validation': {  # 新增
        'N351_hbond': {...},
        'H357_hbond': {...},
        'W400_hbond': {...},
        'total_hbonds': int,
        'is_canonical_gloop': bool
    }
}
```

#### 改进的 `key_residues_engaged` 输出

```python
[
    {
        'position': 'G (pos 6)',
        'residue': 'GLY 65',
        'contacts_glue_vdw': True,
        'min_distance': 3.8,
        'note': '保守Gly面向MGD'
    },
    {
        'position': 'G-4',
        'residue': 'ILE 61',
        'contacts_crbn_sidechain': True,
        'note': '侧链与CRBN接触'
    },
    # ...
]
```

---

## 📊 使用示例

### 示例 1：验证 GSPT1 的 G-loop (PDB 6H0G)

```python
# 在 PyMOL 中
fetch 6H0G, async_=0

# 1. 检测 G-motif
hits = find_crbn_g_motif(
    obj_name='6H0G',
    template_mode='builtin',
    template_builtin='GSPT1 (6H0G A:60-67)',
    rmsd_cutoff=2.0,
    require_gly_pos6=True
)

# 2. 验证 CRBN 氢键
result = validate_crbn_hbonds(
    obj_name='6H0G',
    g_motif_chain='A',
    g_motif_resi_range='60-67',
    crbn_chain='E'
)

# 3. 完整分析（包含 Glue）
binding_result = analyze_g_motif_glue_binding(
    obj_name='6H0G',
    g_motif_chain='A',
    g_motif_resi_range='60-67',
    glue_resname='CC9',  # Lenalidomide 衍生物
    crbn_chain='E',
    validate_hbonds=True
)

print(f"是否为标准 G-loop: {result['is_canonical_gloop']}")
print(f"CRBN 氢键数: {result['total_hbonds']}/3")
print(f"是否为分子胶底物: {binding_result['is_glue_substrate']}")
```

### 示例 2：批量验证候选 G-loop

```python
# 检测所有可能的 G-loop
hits = find_crbn_g_motif(
    obj_name='my_protein',
    rmsd_cutoff=3.5,
    require_gly_pos6=True
)

print(f"发现 {len(hits)} 个候选 G-loop")

# 逐个验证氢键
for i, (chain, start, end, seq, rmsd) in enumerate(hits):
    print(f"\n候选 #{i+1}: {chain}:{start}-{end} {seq}")
    
    hb = validate_crbn_hbonds(
        obj_name='my_protein',
        g_motif_chain=chain,
        g_motif_resi_range=f"{start}-{end}",
        crbn_chain='CRBN_CHAIN'
    )
    
    if hb and hb['is_canonical_gloop']:
        print(f"  ✅ 标准 G-loop: {hb['total_hbonds']}/3 氢键")
    else:
        print(f"  ❌ 非标准 G-loop")
```

### 示例 3：GUI 集成

在 `unified_gui.py` 的 G-motif 分析 Tab 中调用：

```python
def run_gmotif_pocket_analysis(self):
    # ... 获取参数 ...
    
    # 检测 G-motif
    hits = find_crbn_g_motif(obj_name, ...)
    
    # 对每个 hit 验证氢键
    for hit in hits:
        hb_result = validate_crbn_hbonds(
            obj_name, hit['chain'], hit['range'], crbn_chain
        )
        
        # 显示结果
        if hb_result['is_canonical_gloop']:
            self.log(f"✅ {hit['chain']}:{hit['range']} - 标准 G-loop")
        else:
            self.log(f"⚠️ {hit['chain']}:{hit['range']} - 氢键不足")
```

---

## 🎯 判断标准

### 标准 G-loop（Canonical G-loop）

满足以下条件之一：
1. **至少 2/3 CRBN 氢键**（N351, H357, W400）
2. RMSD ≤ 2.5 Å（相对真实模板）
3. 保守 Gly 与 MGD 有 vdW 接触（< 4.5 Å）

### 高置信度分子胶底物

同时满足：
1. ✅ 标准 G-loop（氢键验证通过）
2. ✅ Glue 与 G-loop 接触数 ≥ 2
3. ✅ Glue 与 CRBN 接触数 ≥ 2
4. ✅ 保守 Gly 与 MGD vdW 接触

---

## 📈 与现有功能对比

| 功能 | 原版本 | 改进版本 |
|------|--------|----------|
| G-loop 检测 | RMSD 匹配 | ✅ 保留 |
| Gly6 检测 | 简单存在检查 | ✅ vdW 接触验证 |
| CRBN 氢键 | ❌ 无 | ✅ **新增** |
| 侧链接触 | ❌ 无 | ✅ **新增** |
| 文献标准 | 无 | ✅ **Schrödinger 标准** |

---

## 🧪 验证案例

### 已知 G-loop 结构（文献验证）

| 蛋白 | PDB | 链 | 范围 | G-loop 序列 | 预期氢键 |
|------|-----|-----|------|-------------|----------|
| GSPT1 | 6H0G | A | 60-67 | IKIIQGAK | 3/3 |
| CK1α | 5FQD | A | 36-43 | INTGSEQR | 3/3 |
| IKZF1 | 5HXB | B | 144-151 | QALICDGP | 2/3 |
| IKZF3 | 6H0F | A | 143-150 | QAVQGDGP | 3/3 |

### 测试建议

```bash
# 下载测试结构
fetch 6H0G 5FQD 5HXB 6H0F, async_=0

# 运行验证脚本
run /path/to/test_gloop_validation.py
```

---

## 🔧 参数调优

### 氢键距离阈值

- **默认**: 3.5 Å（兼容 Schrödinger 和宽松标准）
- **严格**: 2.8 Å（仅 Schrödinger 标准）
- **宽松**: 4.0 Å（包含弱氢键）

```python
validate_crbn_hbonds(..., max_hbond_dist=2.8)  # 严格模式
```

### vdW 接触阈值

- **Gly-MGD**: 4.5 Å（文献推荐）
- **侧链-CRBN**: 5.0 Å（默认 distance_threshold）

---

## 📖 参考文献

1. **Schrödinger 标准**: Schrödinger Maestro 氢键定义（≤2.8Å，角度>120°）
2. **G-loop 定义**: Petzold, G. et al. (2016) *Nature* 532, 127–130
3. **CRBN 关键残基**: Hanzl, A. et al. (2023) *Cell Rep.* - 深度突变扫描
4. **Annual Review**: "CRBN-MGD Complexes Recognize a Structural Motif for Degradation" (2023)

---

## 💡 未来改进方向

1. **角度验证**: 添加 D-H-A 角度计算（目前仅距离）
2. **MGD 特异性**: 区分不同 IMiD 类型（Thalidomide/Lenalidomide/Pomalidomide）
3. **Neosubstrate 预测**: 计算无药物时的亲和力（ΔΔRG）
4. **机器学习**: 训练分类器预测 degradability

---

**最后更新**: 2025-11-11  
**维护者**: Vesper  
**基于文献**: Annual Review of Pharmacology and Toxicology Vol. 63:2023
