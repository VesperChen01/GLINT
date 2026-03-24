# FoldX 安装指南

GLINT 使用 FoldX 进行蛋白质突变能量计算（ΔΔG 分析）。FoldX 需要用户自行下载并放置到此目录。

## 下载 FoldX

1. 访问 FoldX 官方网站：https://foldxsuite.crg.eu/
2. 注册账号（免费用于学术研究）
3. 登录后下载适合您操作系统的 FoldX 可执行文件

## 安装步骤

### macOS / Linux

1. 下载 FoldX 可执行文件（通常命名为 `foldx_20270131` 或类似）
2. 将可执行文件复制到此目录：
   ```bash
   cp /path/to/downloaded/foldx_20270131 external/foldx/
   ```
3. 添加执行权限：
   ```bash
   chmod +x external/foldx/foldx_20270131
   ```

### Windows

1. 下载 FoldX 可执行文件（通常命名为 `foldx_1_20270131.exe`）
2. 将 `.exe` 文件复制到此目录

## 验证安装

安装完成后，在 GLINT 中运行突变分析功能时，软件会自动检测此目录中的 FoldX 可执行文件。

## 检测优先级

GLINT 按以下顺序检测 FoldX：

1. **项目目录**（推荐）：`external/foldx/` 目录中的可执行文件
2. **环境变量**：`FOLDX` 环境变量指定的路径
3. **系统 PATH**：系统 PATH 中的 `foldx` 命令

## 支持的可执行文件名

- macOS/Linux: `foldx`, `foldx_20270131`, `foldx5`, `foldx_Mac_0`
- Windows: `foldx.exe`, `foldx_1_20270131.exe`, `foldx5.exe`

## 常见问题

**Q: 为什么不直接包含 FoldX？**  
A: FoldX 需要用户注册并同意许可协议，因此无法直接分发。

**Q: 如何确认 FoldX 已正确安装？**  
A: 在终端运行以下命令测试：
```bash
./external/foldx/foldx_20270131 --version
```

**Q: 可以使用系统全局安装的 FoldX 吗？**  
A: 可以。如果您已将 FoldX 添加到系统 PATH 或设置了 FOLDX 环境变量，GLINT 会自动检测。

## 技术支持

如遇到问题，请检查：
- FoldX 文件是否有执行权限（macOS/Linux）
- 文件名是否符合上述支持的命名规范
- FoldX 版本是否兼容（推荐 FoldX 5.0+）

