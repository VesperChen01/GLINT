#!/bin/bash
# macOS App Bundle Build Script for GLINT Installer

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Building GLINT macOS Installer App Bundle..."

# 清理旧的构建
rm -rf "../GLINT Installer.app"

# 创建 .app 目录结构
APP_NAME="GLINT Installer.app"
APP_PATH="../${APP_NAME}"
CONTENTS_DIR="${APP_PATH}/Contents"
MACOS_DIR="${CONTENTS_DIR}/MacOS"
RESOURCES_DIR="${CONTENTS_DIR}/Resources"

mkdir -p "$MACOS_DIR"
mkdir -p "$RESOURCES_DIR"

echo "Created app bundle structure at: $APP_PATH"

# 复制 GLINT 源码到 Resources
echo "Copying GLINT source code..."
cp -r "../glint" "$RESOURCES_DIR/"

# 复制图标
if [ -f "../glint/assets/logo.png" ]; then
    cp "../glint/assets/logo.png" "$RESOURCES_DIR/AppIcon.png"
    echo "Copied icon"
fi

# 复制安装脚本到 Resources
cp "GLINT_Installer.py" "$RESOURCES_DIR/"
echo "Copied installer script"

# 创建启动器脚本
cat > "$MACOS_DIR/launcher" << 'EOF'
#!/bin/bash
# GLINT Installer Launcher

CONDA_EXE=""
if command -v conda &> /dev/null; then
    CONDA_EXE=$(command -v conda)
else
    for p in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" "/opt/miniconda3/bin/conda" "/opt/anaconda3/bin/conda" "/usr/local/bin/conda" "/opt/homebrew/bin/conda" "/opt/homebrew/Caskroom/miniconda/base/bin/conda"; do
        if [ -x "$p" ]; then
            CONDA_EXE="$p"
            break
        fi
    done
fi

if [ -z "$CONDA_EXE" ]; then
  osascript -e 'display alert "Conda Not Found" message "Please install Miniconda first."'
  exit 1
fi

eval "$($CONDA_EXE shell.bash hook)"
CONDA_BASE="$(conda info --base 2>/dev/null)"
PYTHON="${CONDA_BASE}/bin/python"

if [ ! -x "$PYTHON" ]; then
  PYTHON=$(command -v python3)
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RESOURCES_DIR="$(dirname "$DIR")/Resources"
INSTALLER_SCRIPT="${RESOURCES_DIR}/GLINT_Installer.py"

"$PYTHON" "$INSTALLER_SCRIPT"
EOF

chmod +x "$MACOS_DIR/launcher"
echo "Created launcher script"

# 创建 Info.plist
cat > "$CONTENTS_DIR/Info.plist" << 'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleExecutable</key>
    <string>launcher</string>
    <key>CFBundleIconFile</key>
    <string>AppIcon</string>
    <key>CFBundleIdentifier</key>
    <string>com.glint.installer</string>
    <key>CFBundleName</key>
    <string>GLINT Installer</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>CFBundleShortVersionString</key>
    <string>0.1.28</string>
    <key>CFBundleVersion</key>
    <string>1</string>
    <key>LSMinimumSystemVersion</key>
    <string>10.13</string>
    <key>NSHighResolutionCapable</key>
    <true/>
</dict>
</plist>
EOF

echo "Created Info.plist"

# 创建 icns 图标 (如果有 sips 工具)
if [ -f "$RESOURCES_DIR/AppIcon.png" ]; then
    ICONSET_DIR="$RESOURCES_DIR/AppIcon.iconset"
    mkdir -p "$ICONSET_DIR"
    
    # 生成不同尺寸的图标
    sips -z 16 16     "$RESOURCES_DIR/AppIcon.png" --out "$ICONSET_DIR/icon_16x16.png" 2>/dev/null
    sips -z 32 32     "$RESOURCES_DIR/AppIcon.png" --out "$ICONSET_DIR/icon_16x16@2x.png" 2>/dev/null
    sips -z 32 32     "$RESOURCES_DIR/AppIcon.png" --out "$ICONSET_DIR/icon_32x32.png" 2>/dev/null
    sips -z 64 64     "$RESOURCES_DIR/AppIcon.png" --out "$ICONSET_DIR/icon_32x32@2x.png" 2>/dev/null
    sips -z 128 128   "$RESOURCES_DIR/AppIcon.png" --out "$ICONSET_DIR/icon_128x128.png" 2>/dev/null
    sips -z 256 256   "$RESOURCES_DIR/AppIcon.png" --out "$ICONSET_DIR/icon_128x128@2x.png" 2>/dev/null
    sips -z 256 256   "$RESOURCES_DIR/AppIcon.png" --out "$ICONSET_DIR/icon_256x256.png" 2>/dev/null
    sips -z 512 512   "$RESOURCES_DIR/AppIcon.png" --out "$ICONSET_DIR/icon_256x256@2x.png" 2>/dev/null
    sips -z 512 512   "$RESOURCES_DIR/AppIcon.png" --out "$ICONSET_DIR/icon_512x512.png" 2>/dev/null
    sips -z 1024 1024 "$RESOURCES_DIR/AppIcon.png" --out "$ICONSET_DIR/icon_512x512@2x.png" 2>/dev/null
    
    # 转换为 icns
    iconutil -c icns "$ICONSET_DIR" -o "$RESOURCES_DIR/AppIcon.icns" 2>/dev/null
    rm -rf "$ICONSET_DIR"
    
    if [ -f "$RESOURCES_DIR/AppIcon.icns" ]; then
        echo "Created AppIcon.icns"
    fi
fi

echo ""
echo "✅ App bundle created successfully!"
echo "Location: $APP_PATH"
echo ""
echo "To create a DMG installer, run:"
echo "  ./create_dmg.sh"

