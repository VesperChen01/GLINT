#!/bin/bash
# release_glint.sh - GLINT Release Script
# Creates DMG installer package
#
# Usage: ./release_glint.sh [version]
# Example: ./release_glint.sh v0.2.3

set -e

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Configuration
APP_NAME="GLINT Installer.app"
CONTENTS="${APP_NAME}/Contents"
MACOS="${CONTENTS}/MacOS"
RESOURCES="${CONTENTS}/Resources"
BUILD_DIR="build_dmg"
GLINT_LIB="glint"

# Get version
if [ -n "$1" ]; then
    VERSION="$1"
else
    VERSION=$(python3 -c "from glint._version import __version__; print(__version__)" 2>/dev/null || echo "v0.2.3")
fi

DMG_NAME="GLINT_Installer_${VERSION}.dmg"

echo ""
echo -e "${BLUE}🚀 GLINT Release Script${NC}"
echo -e "${BLUE}========================${NC}"
echo -e "Version: ${GREEN}${VERSION}${NC}"
echo ""

# ============================================================================
# STEP 1: Clean up old files
# ============================================================================
echo -e "${BLUE}[1/5] Cleaning up old files...${NC}"
rm -rf "${APP_NAME}" "${BUILD_DIR}" "${DMG_NAME}"
echo -e "${GREEN}   ✅ Cleaned${NC}"

# ============================================================================
# STEP 2: Create app structure
# ============================================================================
echo -e "${BLUE}[2/5] Creating GLINT Installer.app...${NC}"

mkdir -p "${MACOS}"
mkdir -p "${RESOURCES}"

# Copy icon
if [ -f "glint/assets/AppIcon.icns" ]; then
    cp "glint/assets/AppIcon.icns" "${RESOURCES}/AppIcon.icns"
    cp "glint/assets/logo.png" "${RESOURCES}/AppIcon.png" 2>/dev/null || true
    echo -e "${GREEN}   ✅ Icon copied${NC}"
fi

# Create Info.plist
cat > "${CONTENTS}/Info.plist" << EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleExecutable</key>
    <string>launcher</string>
    <key>CFBundleIconFile</key>
    <string>AppIcon</string>
    <key>CFBundleIdentifier</key>
    <string>com.vesper.glint.installer</string>
    <key>CFBundleName</key>
    <string>GLINT Installer</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>CFBundleShortVersionString</key>
    <string>${VERSION}</string>
</dict>
</plist>
EOF
echo -e "${GREEN}   ✅ Info.plist created${NC}"

# Create launcher script
cat > "${MACOS}/launcher" << 'LAUNCHER_EOF'
#!/bin/bash
# GLINT Installer Launcher

CONDA_EXE=""
if command -v conda &> /dev/null; then
    CONDA_EXE=$(command -v conda)
else
    for p in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" "/opt/miniconda3/bin/conda" "/opt/anaconda3/bin/conda" "/usr/local/bin/conda" "/opt/homebrew/bin/conda"; do
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
  PYTHON=$(command -v python)
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RESOURCES_DIR="$(dirname "$DIR")/Resources"
INSTALLER_SCRIPT="${RESOURCES_DIR}/GLINT_Installer.py"

"$PYTHON" "$INSTALLER_SCRIPT"
LAUNCHER_EOF

chmod +x "${MACOS}/launcher"
echo -e "${GREEN}   ✅ Launcher created${NC}"

# Copy installer Python script
if [ -f "macos_installer/GLINT_Installer.py" ]; then
    cp "macos_installer/GLINT_Installer.py" "${RESOURCES}/GLINT_Installer.py"
    echo -e "${GREEN}   ✅ Installer script copied${NC}"
else
    echo -e "${RED}   ❌ macos_installer/GLINT_Installer.py not found!${NC}"
    exit 1
fi

# Copy glint library
echo -e "${BLUE}[3/5] Copying GLINT library...${NC}"
cp -R "${GLINT_LIB}" "${RESOURCES}/"
echo -e "${GREEN}   ✅ GLINT library copied${NC}"

# ============================================================================
# STEP 4: Create DMG
# ============================================================================
echo -e "${BLUE}[4/5] Creating DMG...${NC}"

mkdir -p "${BUILD_DIR}"
cp -R "${APP_NAME}" "${BUILD_DIR}/"

# Create DMG
hdiutil create -volname "GLINT Installer" -srcfolder "${BUILD_DIR}" -ov -format UDZO "${DMG_NAME}"
echo -e "${GREEN}   ✅ DMG created: ${DMG_NAME}${NC}"

# ============================================================================
# STEP 5: Cleanup
# ============================================================================
echo -e "${BLUE}[5/5] Cleaning up...${NC}"
rm -rf "${BUILD_DIR}"
echo -e "${GREEN}   ✅ Cleanup complete${NC}"

echo ""
echo -e "${GREEN}🎉 Release complete!${NC}"
echo -e "   DMG: ${BLUE}${DMG_NAME}${NC}"
echo -e "   App: ${BLUE}${APP_NAME}${NC}"
ls -lh "${DMG_NAME}"

