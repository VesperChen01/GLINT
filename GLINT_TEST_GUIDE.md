# GLINT 快速测试指南

**日期**: 2026-01-08 00:16  
**状态**: 重命名完成，等待测试

---

## 🚀 在 PyMOL 中测试 GLINT

### 方法 1：重启 PyMOL（推荐）

1. **完全退出 PyMOL**
2. **重新启动 PyMOL**
3. **自动加载**：启动脚本会自动加载 GLINT
4. **查看输出**：
   ```
   ✅ GLINT dev paths loaded
   ✅ GLINT v1.x.x loaded
   Type 'glint_gui' to launch GUI
   ```
5. **启动 GUI**：
   ```python
   glint_gui
   ```

---

### 方法 2：当前 PyMOL 会话中手动加载

在 PyMOL 命令行执行：

```python
# 1. 添加路径
import sys
sys.path.insert(0, "/Users/vesper/Desktop/git/GlueTK")

# 2. 导入 glint（注意：不是 gluetk）
import glint

# 3. 初始化插件
glint.__init_plugin__()

# 4. 启动 GUI
glint_gui
```

---

## ✅ 测试清单

### 基础功能测试
- [ ] GUI 启动成功
- [ ] 窗口标题显示 "GLINT"
- [ ] 所有标签页正常显示

### 核心功能测试
```python
# 1. G-motif 检测
fetch 6H0G
find_crbn_g_motif 6H0G

# 2. PPI 分析
ppi_analyze 6H0G, ['A'], ['B']

# 3. Neo-epitope 识别
neo_epitope_find 6H0G, ['A'], ['B'], 'LEN'
```

### Bug 修复验证
- [ ] **密集线条问题**：加载复合物，运行蛋白-配体分析，检查是否还有密集粉色线条
- [ ] **2D 图显示**：生成 2D 图，检查线条是否遮挡原子
- [ ] **残基颜色**：检查负电性/正电性残基颜色是否清晰

---

## 🐛 常见问题

### 问题 1：`ModuleNotFoundError: No module named 'gluetk'`
**原因**：项目已重命名为 `glint`  
**解决**：使用 `import glint` 而不是 `import gluetk`

### 问题 2：`ModuleNotFoundError: No module named 'glint'`
**原因**：路径未正确添加  
**解决**：
```python
import sys
sys.path.insert(0, "/Users/vesper/Desktop/git/GlueTK")
import glint
```

### 问题 3：GUI 启动失败
**原因**：Qt 绑定问题  
**解决**：
```bash
conda install -c conda-forge pyqt -y
```

### 问题 4：命令 `gluetk_gui` 不存在
**原因**：已重命名为 `glint_gui`  
**解决**：使用 `glint_gui` 命令

---

## 📊 预期输出

### 成功加载的输出
```
🧬 GLINT - Glue Interface Analyzer v1.x.x
┌────────────────────────────────────────────────┐
│  Quick Start:                                   │
│    • glint_gui            - Launch GUI          │
│    • help(glint_gui)      - Show help          │
│    • Plugins → GLINT      - Menu access        │
└────────────────────────────────────────────────┘
```

### GUI 窗口标题
- 旧版本：`GlueTK`
- 新版本：`GLINT` ✅

---

## 🔍 验证修复效果

### 1. π-阳离子相互作用
**测试步骤**：
```python
fetch 6H0G
analyze_protein_ligand_interactions 6H0G, ligand_name='LEN'
```
**预期**：π-阳离子相互作用数量合理（不再有大量假阳性）

### 2. 3D 可视化密集线条
**测试步骤**：
```python
fetch 6H0G
visualize_protein_ligand_3d 6H0G, ligand_name='LEN'
```
**预期**：不再有密集的粉色线条

### 3. 2D 图线条遮挡
**测试步骤**：
```python
fetch 6H0G
generate_2d_diagram 6H0G, ligand_name='LEN'
```
**预期**：相互作用线在配体原子下方，不遮挡

### 4. 残基颜色
**测试步骤**：生成 2D 图，检查残基气泡颜色  
**预期**：
- 负电性残基：明显的红色 (#EF9A9A)
- 正电性残基：明显的紫色 (#CE93D8)

---

## 📝 测试报告模板

```markdown
## GLINT 测试报告

**测试日期**: 2026-01-08
**测试者**: [你的名字]
**PyMOL 版本**: 3.1.3

### 基础功能
- [ ] GUI 启动: ✅/❌
- [ ] 窗口标题: ✅/❌
- [ ] 所有标签页: ✅/❌

### Bug 修复验证
- [ ] 密集线条问题: ✅/❌
- [ ] 2D 图显示: ✅/❌
- [ ] 残基颜色: ✅/❌

### 核心功能
- [ ] G-motif 检测: ✅/❌
- [ ] PPI 分析: ✅/❌
- [ ] Neo-epitope: ✅/❌

### 问题记录
[记录遇到的任何问题]

### 总体评价
[成功/失败/部分成功]
```

---

## 🎯 下一步

### 如果测试成功 ✅
1. 提交 Git 更改
2. 更新 GitHub 仓库名
3. 生成新 Logo
4. 开始撰写论文

### 如果测试失败 ❌
1. 记录错误信息
2. 检查 `sys.path`
3. 验证文件结构
4. 联系支持

---

**祝测试顺利！** 🚀

