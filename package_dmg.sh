#!/bin/bash
# package_dmg.sh
# 打包 GlueTK Installer.app 为 DMG 文件

set -e

# 配置
APP_NAME="GlueTK Installer"
DMG_NAME="GlueTK_Installer_v0.1.3-beta"
BUILD_DIR="build_dmg"
SOURCE_APP="GlueTK Installer.app"
GLUETK_LIB="gluetk"

echo "📦 Packaging GlueTK into DMG..."

# 1. 清理并创建构建目录
rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"

# 2. 复制 Installer.app
echo "   Copying ${APP_NAME}.app..."
cp -R "${SOURCE_APP}" "${BUILD_DIR}/${APP_NAME}.app"

# 3. 将 gluetk 源码目录复制到 Installer.app 内部
#    这样 App 就是自包含的，不需要外部文件夹
echo "   Bundling gluetk source code..."
mkdir -p "${BUILD_DIR}/${APP_NAME}.app/Contents/Resources/gluetk"
# 排除 __pycache__ 等垃圾文件
rsync -av --exclude='__pycache__' --exclude='*.pyc' --exclude='.DS_Store' \
      "${GLUETK_LIB}/" "${BUILD_DIR}/${APP_NAME}.app/Contents/Resources/gluetk/"

# 4. 创建 DMG
echo "   Creating disk image..."
rm -f "${DMG_NAME}.dmg"

# 使用 hdiutil 创建 dmg
hdiutil create -volname "GlueTK Installer" \
    -srcfolder "${BUILD_DIR}" \
    -ov -format UDZO \
    "${DMG_NAME}.dmg"

echo ""
echo "✅ DMG created: ${DMG_NAME}.dmg"
echo "   You can distribute this file to other users."
