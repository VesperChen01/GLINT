# GLINT macOS 分发指南

## 两种DMG类型

我们现在提供两种macOS安装方式，满足不同用户需求：

### 1. 安装器DMG（首次安装）
**文件名**: `GLINT_Installer_vv0.1.26-beta.dmg`  
**大小**: 4.9 MB  
**用途**: 首次安装GLINT

#### 特点：
- ✅ 完整的安装向导
- ✅ 自动检测和安装conda环境
- ✅ 安装所有依赖（PyMOL、RDKit等）
- ✅ 配置GLINT到PyMOL
- ✅ 创建应用程序快捷方式

#### 使用方法：
1. 下载并打开 `GLINT_Installer_vv0.1.26-beta.dmg`
2. 双击 "GLINT Installer.app"
3. 按照安装向导完成安装
4. 安装完成后，在 ~/Applications 中找到 GLINT.app

#### 适用人群：
- 首次安装GLINT的用户
- 需要安装依赖环境的用户
- 不熟悉命令行的用户

---

### 2. 拖拽式DMG（快速安装/更新）
**文件名**: `GLINT_vv0.1.26-beta.dmg`  
**大小**: 4.3 MB  
**用途**: 快速安装或更新GLINT

#### 特点：
- ✅ 标准的macOS拖拽安装方式
- ✅ 打开DMG后看到GLINT.app和Applications文件夹
- ✅ 直接拖拽GLINT.app到Applications即可
- ✅ 适合已安装依赖的用户更新版本

#### 使用方法：
1. 下载并打开 `GLINT_vv0.1.26-beta.dmg`
2. 将 GLINT.app 拖拽到 Applications 文件夹
3. 从Applications启动GLINT

#### 注意事项：
⚠️ **需要先运行过安装器DMG**，因为：
- 需要conda环境已安装
- 需要glint环境已创建
- 需要PyMOL等依赖已安装

如果直接使用拖拽式DMG而没有安装依赖，启动时会提示：
```
GLINT Not Installed
Please run the GLINT Installer first to set up the environment and dependencies.
```

#### 适用人群：
- 已经使用过安装器DMG的用户
- 需要更新GLINT版本的用户
- 熟悉macOS应用安装的用户

---

## 分发策略

### 推荐方案：

#### 对于新用户：
1. 提供 **安装器DMG** 作为主要下载选项
2. 在文档中说明这是完整安装包
3. 强调会自动处理所有依赖

#### 对于现有用户（更新）：
1. 提供 **拖拽式DMG** 作为更新选项
2. 说明这是轻量级更新包
3. 提醒需要先安装过完整版

### 下载页面示例：

```markdown
## 下载 GLINT for macOS

### 首次安装（推荐）
📦 [GLINT_Installer_vv0.1.26-beta.dmg](link) (4.9 MB)
- 完整安装包，包含安装向导
- 自动安装所有依赖
- 适合首次使用的用户

### 快速安装/更新
🚀 [GLINT_vv0.1.26-beta.dmg](link) (4.3 MB)
- 拖拽式安装
- 适合已安装过GLINT的用户
- 快速更新到最新版本
```

---

## 技术细节

### 安装器DMG内容：
```
GLINT Installer.app/
├── Contents/
│   ├── MacOS/launcher          # 启动脚本
│   ├── Resources/
│   │   ├── GLINT_Installer.py  # 安装向导GUI
│   │   ├── glint/              # GLINT源码
│   │   └── AppIcon.icns        # 图标
│   └── Info.plist
```

### 拖拽式DMG内容：
```
/Volumes/GLINT/
├── GLINT.app                   # 可拖拽的应用
│   ├── Contents/
│   │   ├── MacOS/GLINT         # 启动脚本
│   │   ├── Resources/
│   │   │   ├── glint/          # GLINT源码
│   │   │   └── AppIcon.icns    # 图标
│   │   └── Info.plist
└── Applications -> /Applications  # 快捷方式
```

### GLINT.app启动流程：
1. 检查conda是否安装
2. 检查glint环境是否存在
3. 激活glint环境
4. 启动PyMOL并加载GLINT
5. 如果环境不存在，提示用户运行安装器

---

## 构建说明

### 构建两种DMG：
```bash
cd macos_installer
./build_all.sh
```

这会创建：
- `GLINT_Installer_vv0.1.26-beta.dmg` - 安装器
- `GLINT_vv0.1.26-beta.dmg` - 拖拽式

### 单独构建：

**只构建安装器DMG：**
```bash
cd macos_installer
./build_app.sh
./create_dmg.sh
```

**只构建拖拽式DMG：**
```bash
cd macos_installer
./build_glint_app.sh
./create_drag_drop_dmg.sh
```

---

## 用户反馈

### 常见问题：

**Q: 我应该下载哪个DMG？**  
A: 如果是第一次安装，下载安装器DMG。如果已经安装过，下载拖拽式DMG更新。

**Q: 拖拽式DMG提示"GLINT Not Installed"？**  
A: 需要先运行安装器DMG来安装依赖环境。

**Q: 两个DMG有什么区别？**  
A: 安装器DMG包含完整的安装向导和依赖安装；拖拽式DMG只包含应用本身，需要依赖已安装。

**Q: 可以只用拖拽式DMG吗？**  
A: 不可以，必须先运行过安装器DMG至少一次，以安装conda环境和依赖。

---

## 总结

✅ **安装器DMG** = 完整安装包（首次使用）  
✅ **拖拽式DMG** = 快速更新包（已有环境）

这种双DMG策略提供了：
- 新用户的简单安装体验
- 现有用户的快速更新方式
- 符合macOS用户习惯的安装方式

