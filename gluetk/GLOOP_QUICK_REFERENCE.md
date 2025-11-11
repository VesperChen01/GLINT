# G-loop 验证快速参考
## Quick Reference Card for CRBN G-loop Validation

---

## 🚀 快速开始 (Quick Start)

### 在 PyMOL 中使用 (In PyMOL)

```python
# 1. 加载插件
run /path/to/gluetk/__init__.py

# 2. 快速测试 - GSPT1 经典案例
fetch 6H0G, async_=0
validate_crbn_hbonds('6H0G', 'A', '60-67', 'E')
```

---

## 📋 三个核心命令 (Three Core Commands)

### 1️⃣ `find_crbn_g_motif()` - 检测 G-loop

**用途**: 在蛋白结构中寻找潜在的 G-loop/G-motif

**参数**:
```python
find_crbn_g_motif(
    obj_name='protein',           # PyMOL 对象名
    template_mode='ideal',        # 'ideal' | 'builtin' | 'selection'
    rmsd_cutoff=3.5,              # RMSD 阈值（Å）
    require_gly_pos6=True,        # 第6位必须是Gly
    out_csv='output.csv',         # 输出CSV
    auto_highlight=1              # 自动高亮
)
```

**返回**: `[(chain, start, end, seq, rmsd), ...]`

**示例**:
```python
hits = find_crbn_g_motif('6H0G', rmsd_cutoff=2.0)
# 输出: [('A', '60', '67', 'IKIIQGAK', 1.23), ...]
```

---

### 2️⃣ `validate_crbn_hbonds()` - 验证氢键 ⭐ NEW

**用途**: 验证 G-loop 与 CRBN 的 3 个关键氢键

**参数**:
```python
validate_crbn_hbonds(
    obj_name='6H0G',
    g_motif_chain='A',
    g_motif_resi_range='60-67',   # 或 (60, 67)
    crbn_chain='E',
    max_hbond_dist=3.5            # 氢键距离阈值
)
```

**返回**:
```python
{
    'N351_hbond': {'found': True, 'distance': 2.75, 'g_pos': -3},
    'H357_hbond': {'found': True, 'distance': 2.82, 'g_pos': -2},
    'W400_hbond': {'found': True, 'distance': 2.68, 'g_pos': -1},
    'total_hbonds': 3,
    'is_canonical_gloop': True
}
```

**判断标准**:
- ✅ `total_hbonds >= 2` → 标准 G-loop
- ⚠️ `total_hbonds < 2` → 非标准

---

### 3️⃣ `analyze_g_motif_glue_binding()` - 完整分析

**用途**: 综合分析 G-loop、CRBN、Glue 三方相互作用

**参数**:
```python
analyze_g_motif_glue_binding(
    obj_name='6H0G',
    g_motif_chain='A',
    g_motif_resi_range='60-67',
    glue_resname='CC9',           # 分子胶名称
    crbn_chain='E',
    distance_threshold=5.0,       # 接触距离阈值
    validate_hbonds=True          # 启用氢键验证 ⭐
)
```

**返回**:
```python
{
    'is_glue_substrate': True,
    'binding_mode': 'glue-induced',  # 或 'direct' / 'partial'
    'confidence': 0.85,
    'glue_contacts_g_motif': 5,
    'glue_contacts_crbn': 6,
    'crbn_hbond_validation': {...},   # validate_crbn_hbonds 的结果
    'key_residues_engaged': [...]
}
```

---

## 🎯 常见工作流 (Common Workflows)

### 工作流 1: 未知结构筛选 G-loop

```python
# Step 1: 加载结构
fetch 1ABC, async_=0

# Step 2: 检测候选 G-loop
hits = find_crbn_g_motif('1ABC', rmsd_cutoff=3.5)

# Step 3: 逐个验证氢键
for chain, start, end, seq, rmsd in hits:
    print(f"\n候选: {chain}:{start}-{end} {seq}")
    result = validate_crbn_hbonds('1ABC', chain, f"{start}-{end}", 'CRBN_CHAIN')
    
    if result and result['is_canonical_gloop']:
        print(f"  ✅ 标准 G-loop ({result['total_hbonds']}/3 氢键)")
    else:
        print(f"  ❌ 非标准")
```

### 工作流 2: 已知 G-loop 深度分析

```python
# 直接验证氢键
hb = validate_crbn_hbonds('6H0G', 'A', '60-67', 'E')

# 检查 Glue 参与
if hb['is_canonical_gloop']:
    glue_result = analyze_g_motif_glue_binding(
        '6H0G', 'A', '60-67', 'CC9', 'E', validate_hbonds=False  # 已验证
    )
    print(f"分子胶底物: {glue_result['is_glue_substrate']}")
```

### 工作流 3: 批量测试（使用测试脚本）

```python
# 在 PyMOL 中
run /path/to/test_gloop_validation.py

# 快速测试 GSPT1
quick_test_6h0g()

# 或运行完整测试套件
run_all_tests()
```

---

## 📊 解读结果 (Interpreting Results)

### 氢键验证

| 氢键数 | 评级 | 说明 |
|--------|------|------|
| 3/3 | ⭐⭐⭐ | 完美匹配，高置信度 G-loop |
| 2/3 | ⭐⭐ | 标准 G-loop，可接受 |
| 1/3 | ⚠️ | 弱结合，可能非标准 |
| 0/3 | ❌ | 不是 G-loop |

### 关键残基（CRBN 侧）

- **Asn351** (G-3 接受者): 突变为 Asp 废除结合
- **His357** (G-2 接受者): 耐药突变热点
- **Trp400** (G-1 接受者): 必需残基

### Glue 结合模式

| 模式 | 说明 | Glue 作用 |
|------|------|-----------|
| `glue-induced` | ✅ 分子胶机制 | 胶连 G-loop 和 CRBN |
| `direct` | ⚠️ 直接结合 | Glue 不参与 |
| `partial` | 🤔 部分接触 | 需进一步验证 |
| `no-binding` | ❌ 无结合 | 不是底物 |

---

## 🔧 参数调优 (Parameter Tuning)

### 氢键距离阈值

```python
# 严格模式（Schrödinger 标准）
validate_crbn_hbonds(..., max_hbond_dist=2.8)

# 默认（兼容）
validate_crbn_hbonds(..., max_hbond_dist=3.5)  # ← 推荐

# 宽松（包含弱氢键）
validate_crbn_hbonds(..., max_hbond_dist=4.0)
```

### RMSD 阈值（G-loop 检测）

```python
# 严格（仅高度相似）
find_crbn_g_motif(..., rmsd_cutoff=2.0)

# 默认
find_crbn_g_motif(..., rmsd_cutoff=3.5)  # ← 推荐

# 宽松（探索性）
find_crbn_g_motif(..., rmsd_cutoff=5.0)
```

---

## 📖 已知案例 (Known Cases)

### 文献验证的 G-loop

| 蛋白 | PDB | G-loop | CRBN链 | Glue | 氢键 | 说明 |
|------|-----|--------|--------|------|------|------|
| GSPT1 | 6H0G | A:60-67 | E | CC9 | 3/3 | 经典案例 |
| CK1α | 5FQD | A:36-43 | B | LEN | 3/3 | Lenalidomide |
| IKZF1 | 5HXB | B:144-151 | C | POM | 2/3 | Pomalidomide |
| IKZF3 | 6H0F | A:143-150 | E | CC9 | 3/3 | Aiolos |

### 测试命令

```python
# GSPT1 - 经典案例
fetch 6H0G, async_=0
validate_crbn_hbonds('6H0G', 'A', '60-67', 'E')

# CK1α - Lenalidomide
fetch 5FQD, async_=0
validate_crbn_hbonds('5FQD', 'A', '36-43', 'B')

# IKZF3 - Aiolos
fetch 6H0F, async_=0
validate_crbn_hbonds('6H0F', 'A', '143-150', 'E')
```

---

## ⚠️ 常见问题 (Troubleshooting)

### 问题 1: 找不到氢键（返回 0/3）

**可能原因**:
- CRBN 链 ID 错误
- G-loop 残基范围不正确
- 阈值过于严格

**解决方案**:
```python
# 检查链 ID
cmd.iterate('obj_name and polymer', 'print(chain)')

# 放宽阈值
validate_crbn_hbonds(..., max_hbond_dist=4.0)
```

### 问题 2: G-loop 检测无结果

**可能原因**:
- Gly6 要求过严（某些 G-loop 不含 Gly）
- RMSD 阈值过低

**解决方案**:
```python
# 关闭 Gly 要求
find_crbn_g_motif(..., require_gly_pos6=False)

# 放宽 RMSD
find_crbn_g_motif(..., rmsd_cutoff=5.0)
```

### 问题 3: 氢键验证返回 None

**可能原因**:
- G-loop 少于 8 个残基
- CRBN 关键残基缺失

**解决方案**:
```python
# 检查残基数
cmd.count_atoms('obj and chain A and resi 60-67')

# 检查 CRBN 关键残基
cmd.iterate('obj and chain E and resi 351+357+400', 'print(resn, resi)')
```

---

## 📚 文档链接 (Documentation)

- **详细指南**: `G_LOOP_VALIDATION_GUIDE.md`
- **测试脚本**: `test_gloop_validation.py`
- **口袋分析**: `GUI_ADVANCED_POCKET_GUIDE.md`
- **项目文档**: `WARP.md`

---

## 💡 专业提示 (Pro Tips)

1. **先检测，后验证**: 用 `find_crbn_g_motif()` 找候选，再用 `validate_crbn_hbonds()` 确认
2. **置信度优先**: ≥2/3 氢键的 G-loop 才是可靠的
3. **Glue 必须**: 真正的分子胶底物必须有 Glue 参与（`glue-induced` 模式）
4. **保守 Gly**: 第 6 位 Gly 与药物的 vdW 接触是关键特征
5. **对比文献**: 用已知案例（GSPT1/6H0G）校准你的参数

---

**版本**: 1.0.0  
**更新**: 2025-11-11  
**文献**: Annual Review of Pharmacology and Toxicology Vol. 63:2023
