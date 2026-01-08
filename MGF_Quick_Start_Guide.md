# 🚀 Molecular Glue Fingerprint (MGF) 快速开始指南

## 📋 目录
- [安装](#安装)
- [5分钟快速入门](#5分钟快速入门)
- [基础用法](#基础用法)
- [高级功能](#高级功能)
- [常见问题](#常见问题)

---

## 🔧 安装

### 前置要求
```bash
# Python 3.8+
# PyMOL 2.5+
# GLINT平台已安装
```

### 安装MGF模块
```bash
cd /path/to/GLINT
git pull origin main  # 获取最新代码

# 安装依赖
pip install numpy pandas matplotlib seaborn scikit-learn scipy
```

### 在PyMOL中加载
```python
# 在PyMOL命令行中
run /path/to/GLINT/glint/molecular_glue_fingerprint.py
```

---

## ⚡ 5分钟快速入门

### 示例：分析Lenalidomide-CRBN-CK1α复合物

```python
# 1. 加载结构
fetch 5fqd

# 2. 生成分子胶指纹
generate_molecular_glue_fingerprint 5fqd, A+B, C, 1NH

# 3. 3D可视化
visualize_molecular_glue_fingerprint_3d 5fqd, result, 1NH
```

**就这么简单！** 🎉

---

## 📖 基础用法

### 1. 生成指纹

```python
# 语法
generate_molecular_glue_fingerprint <object>, <e3_chains>, <substrate_chains>, <glue_resname>

# 示例1: 单链E3，单链底物
generate_molecular_glue_fingerprint 6h0g, A, B, CC9

# 示例2: 多链E3（CRBN是二聚体）
generate_molecular_glue_fingerprint 5fqd, A+B, C, 1NH

# 示例3: 保存结果到文件
generate_molecular_glue_fingerprint 6h0g, A, B, CC9, output.json
```

### 2. 解读输出

```
============================================================
  Molecular Glue Fingerprint Analysis
============================================================
E3 Chains: ['A']
Substrate Chains: ['B']
Glue: CC9

--- Interface Statistics ---
E3-Glue interactions: 15          ← E3与分子胶的相互作用数
Glue-Substrate interactions: 12   ← 分子胶与底物的相互作用数
Neo-epitope contacts: 18          ← E3与底物的直接接触数

Bridging atoms: 5                 ← 桥接原子数（越多越好）
Cooperativity score: 0.724        ← 协同性评分（0-1，越高越好）

--- Bridging Atoms (Key for Glue Function) ---
  1. C10 (resi 1501)
  2. C14 (resi 1501)
  3. N3 (resi 1501)
  4. O2 (resi 1501)
  5. C7 (resi 1501)
```

### 3. 评分解读

| 指标 | 优秀 | 良好 | 一般 | 差 |
|-----|------|------|------|-----|
| **协同性评分** | >0.7 | 0.5-0.7 | 0.3-0.5 | <0.3 |
| **桥接原子数** | >5 | 3-5 | 1-3 | 0 |
| **Neo-epitope接触** | >15 | 10-15 | 5-10 | <5 |
| **总相互作用数** | >25 | 15-25 | 10-15 | <10 |

---

## 🎨 可视化

### 1. 3D可视化（PyMOL）

```python
# 自动高亮关键特征
visualize_molecular_glue_fingerprint_3d 5fqd, result, 1NH

# 图例：
# 🟡 黄色球体 = 桥接原子（关键！）
# 🔵 青色棒状 = E3-Glue界面残基
# 🔴 洋红棒状 = Glue-底物界面残基
# 🟠 橙色表面 = Neo-epitope（新界面）
```

### 2. 热图（Python脚本）

```python
from glint.molecular_glue_fingerprint import plot_molecular_glue_fingerprint_heatmap
import json

# 加载指纹数据
with open('output.json', 'r') as f:
    fp_result = json.load(f)

# 生成热图
plot_molecular_glue_fingerprint_heatmap(fp_result, 'heatmap.png')
```

### 3. 雷达图比较多个化合物

```python
from glint.molecular_glue_fingerprint import plot_cooperativity_radar

# 加载多个指纹
fp_results = [fp1, fp2, fp3]
compound_names = ['Compound A', 'Compound B', 'Compound C']

plot_cooperativity_radar(fp_results, compound_names, 'radar.png')
```

---

## 🔬 高级功能

### 1. 批量分析

```python
# batch_analysis.py
from pymol import cmd
from glint.molecular_glue_fingerprint import MolecularGlueFingerprint
import pandas as pd

pdb_files = ['compound1.pdb', 'compound2.pdb', 'compound3.pdb']
results = []

for pdb in pdb_files:
    obj_name = pdb.replace('.pdb', '')
    cmd.load(pdb, obj_name)
    
    mgf = MolecularGlueFingerprint(obj_name, ['A'], ['B'], 'LIG')
    fp = mgf.generate_fingerprint()
    
    results.append({
        'compound': obj_name,
        'cooperativity': fp['cooperativity_score'],
        'bridging_atoms': len(fp['bridging_atoms']),
        'e3_glue_interactions': fp['metadata']['n_e3_glue_interactions'],
        'glue_sub_interactions': fp['metadata']['n_glue_substrate_interactions']
    })
    
    cmd.delete(obj_name)

df = pd.DataFrame(results)
df.to_csv('batch_results.csv', index=False)
print(df.sort_values('cooperativity', ascending=False))
```

### 2. SAR分析

```python
from glint.molecular_glue_fingerprint import compare_molecular_glue_analogs

# 比较类似物
sar_results = compare_molecular_glue_analogs(
    reference_pdb='lead_compound.pdb',
    analog_pdbs=['analog1.pdb', 'analog2.pdb', 'analog3.pdb'],
    e3_chains=['A'],
    substrate_chains=['B'],
    glue_resname='LIG'
)

# 查看SAR洞察
for insight in sar_results['sar_insights']:
    print(insight)
```

### 3. MD轨迹分析

```python
from glint.molecular_glue_fingerprint_md import MolecularGlueFingerprintMD

# 分析MD轨迹中的指纹演化
mgf_md = MolecularGlueFingerprintMD(
    e3_chains=['A'],
    substrate_chains=['B'],
    glue_resname='LIG'
)

md_result = mgf_md.run_trajectory(
    trajectory_file='md_trajectory.dcd',
    topology_file='topology.pdb',
    frame_step=10  # 每10帧采样一次
)

# 可视化时间演化
from glint.molecular_glue_fingerprint import plot_fingerprint_evolution
plot_fingerprint_evolution(md_result, 'md_evolution.png')

print(f"平均协同性: {np.mean(md_result['cooperativity_timeline']):.3f}")
print(f"稳定性评分: {md_result['stability_score']:.3f}")
```

### 4. 虚拟筛选排序

```python
from glint.molecular_glue_fingerprint import rank_molecular_glue_candidates

# 对候选化合物排序
pdb_files = [f'candidate_{i}.pdb' for i in range(1, 101)]

ranked_df = rank_molecular_glue_candidates(
    pdb_files=pdb_files,
    e3_chains=['A'],
    substrate_chains=['B'],
    glue_resname='LIG'
)

# 查看Top 10
print(ranked_df.head(10))

# 保存结果
ranked_df.to_csv('virtual_screening_results.csv', index=False)
```

---

## 💡 实用技巧

### 技巧1: 快速识别链ID

```python
# 在PyMOL中查看链信息
fetch 5fqd
print(cmd.get_chains('5fqd'))  # 输出: ['A', 'B', 'C']

# 查看配体残基名
iterate 5fqd and organic, print(resn)
```

### 技巧2: 处理多个配体

```python
# 如果结构中有多个配体分子
# 使用resi（残基编号）指定特定配体

# 方法1: 先选择特定配体
select my_glue, 5fqd and resn 1NH and resi 1501

# 方法2: 在分析时指定
# 修改代码以支持resi参数
```

### 技巧3: 导出高质量图片

```python
# PyMOL中导出高分辨率图片
ray 2400, 2400  # 设置分辨率
png molecular_glue_3d.png, dpi=300

# Python脚本中设置
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
```

### 技巧4: 自定义协同性计算

```python
# 如果需要自定义协同性评分公式
class CustomMGF(MolecularGlueFingerprint):
    def _calculate_cooperativity(self):
        # 自定义公式
        e3_glue_strength = self._e3_glue_interactions.total_count()
        glue_sub_strength = self._glue_substrate_interactions.total_count()
        neo_strength = len(self._neo_epitope_interactions.get('interface_residues', []))

        # 例如：强调Neo-epitope的重要性
        cooperativity = (neo_strength * 2) / (e3_glue_strength + glue_sub_strength + neo_strength)

        return cooperativity
```

---

## ❓ 常见问题

### Q1: 如何确定E3和底物的链ID？

**A**: 使用以下方法：

```python
# 方法1: PyMOL GUI
# 点击结构 → 右侧显示链信息

# 方法2: 命令行
fetch 5fqd
iterate 5fqd and name CA, print(f"{chain}: {resn}{resi}")

# 方法3: 查看PDB文件头部
# COMPND行通常包含链信息
```

### Q2: 协同性评分为0是什么原因？

**A**: 可能的原因：
1. **Neo-epitope不存在**：E3和底物没有直接接触
2. **链ID错误**：检查E3和底物链是否正确
3. **距离阈值过严**：尝试增大cutoff参数
4. **结构质量问题**：检查PDB文件完整性

```python
# 调试方法
mgf = MolecularGlueFingerprint('5fqd', ['A', 'B'], ['C'], '1NH')
result = mgf.generate_fingerprint()

print(f"E3-Glue: {result['metadata']['n_e3_glue_interactions']}")
print(f"Glue-Sub: {result['metadata']['n_glue_substrate_interactions']}")
print(f"Neo-epitope: {result['metadata']['n_neo_epitope_contacts']}")
```

### Q3: 如何处理缺失的原子或残基？

**A**:
```python
# 1. 使用PyMOL修复
remove 5fqd and not polymer  # 移除非聚合物
h_add 5fqd  # 添加氢原子

# 2. 使用MODELLER补全缺失残基
# 3. 使用PDBFixer
from pdbfixer import PDBFixer
fixer = PDBFixer(filename='input.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)
```

### Q4: 桥接原子数量很少怎么办？

**A**: 桥接原子少可能表明：
1. **分子胶效率低**：这是真实的生物学信号
2. **距离阈值需调整**：尝试增大cutoff
3. **结构构象不佳**：考虑MD模拟优化

```python
# 调整距离阈值
mgf._is_in_contact(coord, chains, cutoff=5.0)  # 默认4.5Å
```

### Q5: 如何比较不同类型的分子胶？

**A**: 使用归一化指纹：

```python
from sklearn.preprocessing import StandardScaler

# 收集多个分子胶的指纹
fingerprints = [fp1['combined_fp'], fp2['combined_fp'], fp3['combined_fp']]

# 标准化
scaler = StandardScaler()
normalized_fps = scaler.fit_transform(fingerprints)

# 计算相似度
from sklearn.metrics.pairwise import cosine_similarity
similarity_matrix = cosine_similarity(normalized_fps)
```

### Q6: 能否用于PROTAC分析？

**A**: 可以！PROTAC也是三元复合物：

```python
# PROTAC = E3-PROTAC-POI
generate_molecular_glue_fingerprint protac_complex, E3_chains, POI_chains, PROTAC_resname

# 注意：PROTAC通常有更长的linker，桥接原子可能更多
```

---

## 📚 示例数据集

### 已验证的分子胶结构

| PDB ID | 分子胶 | E3 | 底物 | 链配置 | 命令示例 |
|--------|--------|----|----|--------|---------|
| 5FQD | Lenalidomide | CRBN | CK1α | E3:A,B Sub:C | `generate_molecular_glue_fingerprint 5fqd, A+B, C, 1NH` |
| 6H0G | CC-885 | CRBN | GSPT1 | E3:A Sub:B | `generate_molecular_glue_fingerprint 6h0g, A, B, CC9` |
| 7S4K | Indisulam | DCAF15 | RBM39 | E3:A,B Sub:C | `generate_molecular_glue_fingerprint 7s4k, A+B, C, 9KP` |
| 7JTO | E7820 | DCAF15 | RBM23 | E3:A,B Sub:C | `generate_molecular_glue_fingerprint 7jto, A+B, C, LXE` |

### 下载示例数据

```bash
# 下载所有示例结构
pymol -c -d "fetch 5fqd 6h0g 7s4k 7jto; save examples.pse"
```

---

## 🎯 最佳实践

### 1. 工作流程建议

```
1. 结构准备
   ↓
2. 链ID确认
   ↓
3. 生成指纹
   ↓
4. 3D可视化验证
   ↓
5. 导出数据分析
   ↓
6. 批量比较（如需要）
```

### 2. 质量控制检查清单

- [ ] 结构完整性：无缺失残基
- [ ] 配体位置：在结合口袋中
- [ ] 链ID正确：E3和底物不混淆
- [ ] 相互作用合理：氢键距离<3.5Å
- [ ] 协同性评分：>0.3为有效分子胶

### 3. 性能优化

```python
# 对于大规模筛选
import multiprocessing as mp

def analyze_single(pdb_file):
    # 分析单个结构
    pass

# 并行处理
with mp.Pool(processes=8) as pool:
    results = pool.map(analyze_single, pdb_files)
```

---

## 🔗 相关资源

- **GLINT文档**: https://glint.readthedocs.io/
- **ProLIF文档**: https://prolif.readthedocs.io/
- **PyMOL Wiki**: https://pymolwiki.org/
- **PDB数据库**: https://www.rcsb.org/

---

## 📞 获取帮助

遇到问题？

1. **查看文档**: `ProLIF_GLINT_Integration_Analysis.md`
2. **GitHub Issues**: 提交bug报告
3. **社区讨论**: 加入GLINT用户群
4. **邮件支持**: support@glint-platform.org

---

**祝您使用愉快！** 🚀

*最后更新: 2026-01-08*


