# 批量热图功能测试示例

## 测试数据格式

每个 CSV 文件代表一个受体（蛋白）与多个配体的结合能数据。

### 文件命名规范
`<受体名称>_scores.csv`

例如：
- `CRBN_scores.csv`
- `VHL_scores.csv`
- `MDM2_scores.csv`

### CSV 内容格式

```csv
Ligand,BindingEnergy
Compound_A,-7.5
Compound_B,-6.8
Compound_C,-8.2
```

或者：

```csv
Name,Affinity
LIG1,-7.5
LIG2,-6.8
LIG3,-8.2
```

## 使用方法

### 1. GUI 方式

1. 启动 PyMOL 插件：`molstruct_gui`
2. 切换到 **Scoring** 标签页
3. 滚动到 **🔥 Batch Heatmap** 卡片
4. 点击 📁 浏览按钮，选择包含 CSV 文件的文件夹
5. （可选）修改文件匹配模式（默认 `*_scores.csv`）
6. 点击 **📊 Generate Heatmap** 按钮
7. 等待生成完成，会弹出对话框询问是否打开图片

### 2. 命令行方式

```python
# 从文件夹生成热图
plot_binding_heatmap /path/to/csv_folder

# 指定输出路径
plot_binding_heatmap /path/to/csv_folder, /path/to/output.png
```

## 输出结果

生成的热图将包含：
- **X 轴**：受体（蛋白）名称
- **Y 轴**：配体名称
- **颜色**：蓝色渐变（深蓝 = 强结合）
- **数值**：结合能（kcal/mol）

热图自动保存为 `binding_energy_heatmap.png`（高分辨率 300 DPI）
