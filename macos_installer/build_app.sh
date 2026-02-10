#!/bin/bash
# macOS App Bundle Build Script for GLINT Installer
set -euo pipefail

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

GLINT_SOURCE_DIR="../glint"
if [ ! -d "$GLINT_SOURCE_DIR" ]; then
    echo "❌ Error: GLINT source directory not found: $GLINT_SOURCE_DIR"
    echo "   Please ensure project structure includes ./glint at repository root."
    exit 1
fi

cp -R "$GLINT_SOURCE_DIR" "$RESOURCES_DIR/glint"
echo "Copied GLINT source from: $GLINT_SOURCE_DIR"

# 复制图标
if [ -f "$GLINT_SOURCE_DIR/assets/logo.png" ]; then
    cp "$GLINT_SOURCE_DIR/assets/logo.png" "$RESOURCES_DIR/AppIcon.png"
    echo "Copied icon"
fi

# 复制安装脚本到 Resources
cp "GLINT_Installer.py" "$RESOURCES_DIR/"
echo "Copied installer script"

# 创建启动器脚本
cat > "$MACOS_DIR/launcher" << 'EOF'
#!/bin/bash
# GLINT Installer Launcher

# Finder launches apps with a very limited PATH.
# Normalize PATH first, then choose a Python with tkinter support.
export PATH="/Library/Frameworks/Python.framework/Versions/Current/bin:/usr/local/bin:/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin"

LOG_DIR="$HOME/Library/Logs"
LOG_FILE="$LOG_DIR/GLINT_Installer.log"
mkdir -p "$LOG_DIR"

PYTHON=""
CANDIDATES=(
  "/Library/Frameworks/Python.framework/Versions/Current/bin/python3"
  "/Library/Frameworks/Python.framework/Versions/3.13/bin/python3"
  "/Library/Frameworks/Python.framework/Versions/3.12/bin/python3"
  "$(command -v python3 || true)"
  "/usr/local/bin/python3"
  "/opt/homebrew/bin/python3"
  "/usr/bin/python3"
)

for p in "${CANDIDATES[@]}"; do
  if [ -n "$p" ] && [ -x "$p" ]; then
    if "$p" -c "import tkinter" >/dev/null 2>&1; then
      PYTHON="$p"
      break
    fi
  fi
done

if [ -z "$PYTHON" ]; then
  osascript -e 'display alert "Tkinter Not Available" message "GLINT Installer requires Python with tkinter GUI support.\n\nPlease install Python from python.org (3.12/3.13 recommended), then relaunch installer.\n\nDetails in ~/Library/Logs/GLINT_Installer.log"'
  {
    echo "[$(date)] tkinter not available in detected python interpreters."
    printf 'Checked candidates:\n'; printf '  %s\n' "${CANDIDATES[@]}"
  } >>"$LOG_FILE"
  exit 1
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RESOURCES_DIR="$(dirname "$DIR")/Resources"
INSTALLER_SCRIPT="${RESOURCES_DIR}/GLINT_Installer.py"

{
  echo "[$(date)] Launching installer"
  echo "Python: $PYTHON"
  "$PYTHON" -c 'import sys, tkinter as tk; print("Python:", sys.version); print("Tk:", tk.TkVersion, "Tcl:", tk.TclVersion)' 2>&1
} >>"$LOG_FILE"

"$PYTHON" "$INSTALLER_SCRIPT" >>"$LOG_FILE" 2>&1
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  osascript -e 'display alert "GLINT Installer Error" message "Failed to launch installer UI. See ~/Library/Logs/GLINT_Installer.log for details."'
  exit $EXIT_CODE
fi
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
    <string>0.2.1</string>
    <key>CFBundleVersion</key>
    <string>1</string>
    <key>LSMinimumSystemVersion</key>
    <string>10.13</string>
    <key>NSRequiresAquaSystemAppearance</key>
    <true/>
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

