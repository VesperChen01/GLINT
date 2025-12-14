# HADDOCK3 故障排除指南

## 问题: CNS 库依赖错误

### 错误信息
```
dyld: Library not loaded: /opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.0.dylib
Referenced from: .../haddock/bin/cns
Reason: no such file
```

### 原因
HADDOCK3 使用 CNS (Crystallography & NMR System) 引擎,CNS 依赖于 GCC 编译器的运行时库 `libquadmath.0.dylib`。

### 解决方案

#### 方法 1: 安装 GCC (推荐)

```bash
# 安装 GCC
brew install gcc

# 验证安装
brew list gcc | grep libquadmath

# 如果找不到符号链接,手动创建
# 首先找到实际的 lib 目录
GCC_LIB=$(brew --prefix gcc)/lib/gcc/current

# 如果 current 不存在,找到具体版本
GCC_VERSION=$(ls $(brew --prefix gcc)/lib/gcc/ | head -1)
GCC_LIB=$(brew --prefix gcc)/lib/gcc/${GCC_VERSION}

# 验证库文件存在
ls -la ${GCC_LIB}/libquadmath.0.dylib
```

#### 方法 2: 使用 Conda 版本的 HADDOCK3 (更稳定)

```bash
# 激活 gluetk 环境
conda activate gluetk

# 卸载 pip 版本
pip uninstall haddock3 -y

# 安装 conda 版本 (预编译,更兼容 macOS)
conda install -c conda-forge -c haddocking haddock3 -y
```

#### 方法 3: 重新构建安装脚本

如果你还没有安装,直接使用更新后的安装脚本:

```bash
./install_gluetk.sh
```

新版本的安装脚本会自动:
1. 检查 GCC 是否已安装
2. 如果未安装则自动安装 GCC
3. 优先使用 conda 安装 HADDOCK3

### 验证修复

在 conda 环境中测试 HADDOCK3:

```bash
conda activate gluetk

# 测试 haddock3 命令
haddock3 -v

# 测试 CNS
python -c "from haddock.modules.topology.topoaa import HaddockModule; print('CNS OK')"
```

## macOS 版本兼容性

**如果遇到 "built for macOS 15.0" 错误:**

这表示 HADDOCK3/CNS 是为较新的 macOS 编译的。解决方案:

1. **升级 macOS** (推荐)
   - 升级到 macOS 15 (Sequoia) 或更新版本

2. **使用兼容版本**
   ```bash
   # 安装旧版本 HADDOCK3
   conda install -c conda-forge -c haddocking "haddock3<3.0" -y
   ```

3. **使用 Docker 容器** (高级用户)
   ```bash
   # 使用 HADDOCK3 官方 Docker 镜像
   docker pull haddocking/haddock3:latest
   ```

## 常见问题

### Q: 为什么需要 GCC?
A: HADDOCK3 的 CNS 引擎是用 Fortran 编译的,需要 GCC 的运行时库(`libgfortran`, `libquadmath`)。

### Q: 可以不用 HADDOCK3 吗?
A: 可以。HADDOCK3 仅用于蛋白-蛋白对接功能。GlueTK 的其他功能(分子胶分析、口袋检测、Vina 对接等)不依赖 HADDOCK3。

### Q: 安装 GCC 需要多久?
A: 通常需要 5-10 分钟,取决于网络速度和电脑性能。

### Q: 如何卸载 HADDOCK3?
```bash
conda activate gluetk
# 如果是 pip 安装的
pip uninstall haddock3
# 如果是 conda 安装的  
conda remove haddock3
```

## 联系支持

如果以上方法都无法解决问题:

1. 查看 HADDOCK3 官方文档: https://github.com/haddocking/haddock3
2. 提交 Issue: https://github.com/haddocking/haddock3/issues
3. 或在 GlueTK 仓库报告问题

## 更新日志

- **v0.1.4-beta**: 添加自动 GCC 检查和安装
- 优先使用 conda 安装 HADDOCK3
- 改进错误提示和路径处理
