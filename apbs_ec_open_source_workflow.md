# 用 APBS + 开源脚本实现 Electrostatic Complementarity（EC）分析

> 目标：在**不依赖商业软件**的前提下，复现类似 Flare 的 **PIP + EC map + EC score** 工作流，并可扩展到\*\*分子胶（三元复合物）\*\*体系。

---

## 1. 核心思想（与 Flare 的一一对应）

- **PIP（Protein Interaction Potential）≈ 蛋白单独的 PB 电势场**\
  使用 Poisson–Boltzmann（PB）方程求解蛋白在空间中的静电势： \(\phi_\text{protein}(x)\)

- **EC（Electrostatic Complementarity）**\
  在配体（或分子胶）**表面采样点** \(x_i\) 上，比较：

  - 蛋白给出的电势 \(\phi_\text{protein}(x_i)\)
  - 配体自身的电势 \(\phi_\text{ligand}(x_i)\)

  同号 → 冲突；异号且幅度匹配 → 互补。

---

## 2. 软件与工具（全部开源）

- **PDB2PQR**：蛋白加氢、赋电荷与半径，生成 `.pqr`
- **APBS**：求解 PB 方程，输出电势网格（OpenDX `.dx`）
- **RDKit**：生成配体 3D、电荷（Gasteiger）
- **NanoShaper / 自定义采样**：生成配体表面点
- **Python (NumPy)**：插值、打分、统计
- **PyMOL / ChimeraX**：可视化 EC map

---

## 3. Step-by-Step 实现流程

### Step 0｜结构准备（SOP 必须统一）

- 固定：pH、His 质子化型、Asp/Glu/Lys/Arg 状态
- 决定是否保留关键水分子
- 配体/分子胶：固定电离态与互变异构体

> ⚠️ 如果比较一系列类似物，**所有准备规则必须一致**。

---

### Step 1｜PDB → PQR（PDB2PQR）

```bash
# 仅保留蛋白（去掉配体、离子）
pdb2pqr --ff=AMBER --with-ph=7.4 protein.pdb protein.pqr
```

输出：`protein.pqr`

---

### Step 2｜APBS 计算蛋白电势（PIP 等价物）

示例 `protein.in`：

```text
read
  mol pqr protein.pqr
end

elec
  mg-auto
  mol 1
  lpbe
  bcfl sdh
  pdie 2.0
  sdie 78.54
  chgm spl2
  srfm smol
  swin 0.3
  sdens 10.0
  temp 298.15
  calcenergy no
  calcforce no

  dime 161 161 161
  cglen 60 60 60
  fglen 40 40 40
  cgcent mol 1
  fgcent mol 1

  write pot dx protein_pot
end
quit
```

运行：

```bash
apbs protein.in
```

输出：`protein_pot.dx`

---

**推荐方案**：

- 使用 **NanoShaper** 生成 SES/SAS 三角网格
- 顶点坐标即为采样点

**轻量方案**（近似）：

- 用 RDKit 得到配体原子坐标与 vdW 半径
- 在每个原子球面均匀采样点
- 剔除被其他原子遮挡的点

---

### Step 4｜读取 `.dx` 并插值蛋白电势

- APBS 输出为规则 3D 网格
- 对每个表面点 \(x_i\) 做**三线性插值**

Python 核心逻辑：

```python
phi_protein = trilinear_interpolate(origin, step, grid, surface_xyz)
```

---

#### 方案 A（推荐，足够做 SAR）

- RDKit 生成 **Gasteiger 电荷**
- 库仑势近似：

$$
\phi_\text{ligand}(x_i) = \sum_j \frac{q_j}{r_{ij}}
$$

#### 方案 B（更严格）

- 将配体也转成 PQR
- 用 APBS/DelPhi 计算配体 PB 电势

---

### Step 6｜定义 EC 局部函数与 EC score

一个稳定、常用的互补性定义：

$$
EC_i = -\frac{2\,\phi_p\,\phi_l}{\phi_p^2 + \phi_l^2 + \varepsilon}
$$

性质：

- 异号 → 正值（互补）
- 同号 → 负值（冲突）
- 幅度匹配 → 接近 +1

整体分数：

$$
EC_\text{score} = \frac{1}{N}\sum_i EC_i
$$

Python 示例：

```python
def ec_local(phi_p, phi_l, eps=1e-6):
    return -(2.0 * phi_p * phi_l) / (phi_p**2 + phi_l**2 + eps)

phi_p = np.clip(phi_p, -10.0, 10.0)
phi_l = np.clip(phi_l, -10.0, 10.0)

ec_i = ec_local(phi_p, phi_l)
ec_score = ec_i.mean()
```

---

### Step 7｜生成 EC map（可视化）

- 将表面点写成 **伪 PDB**
- 把 `EC_i` 写入 **B-factor**
- 用 PyMOL/ChimeraX 着色

```python
write_points_pdb("ec_points.pdb", surface_xyz, ec_i * 100)
```

PyMOL：

```pymol
load ec_points.pdb
spectrum b, blue_white_red, ec_points, minimum=-100, maximum=100
show spheres, ec_points
```

---

## 4. 分子胶（三元体系）的扩展

至少计算 **2 + 1 个界面**：

1. **EC(A–glue)**：Protein A 单独跑 APBS → 胶表面 EC
2. **EC(B–glue)**：Protein B 单独跑 APBS → 胶表面 EC
3. **EC(A–B)**（可选但推荐）：胶存在下的 PPI 新界面

> 实务建议：
>
> - 不要只用单一构象
> - 对 top-N 对接姿势或 MD 抽帧统计 EC 分布（均值 / 分位数）

---

## 5. 常见坑位（Checklist）

-

---

## 6. 适用场景总结

- 解释 **电性主导的 SAR**（如 CF₃ vs CH₃）
- 指导分子胶取代基设计
- 快速筛选“电势互补更合理”的构象与修饰方案

---

> 本流程本质是：**PB 电势 + 表面采样 + 互补函数**。\
> Flare 把它产品化；你现在有的是**可复现、可扩展、可审计的开源版本**。

