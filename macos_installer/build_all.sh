#!/bin/bash
# Complete build script for GLINT macOS distribution
# Creates both installer DMG and drag-and-drop DMG
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "GLINT macOS - Complete Build"
echo "=========================================="
echo ""

# Get version
VERSION=""
if [ -f ../glint/_version.py ]; then
    VERSION=$(grep "__version__" ../glint/_version.py | cut -d'"' -f2 || true)
fi
if [ -z "$VERSION" ]; then
    VERSION="0.2.3"
fi
VERSION_CLEAN="${VERSION#v}"

echo "Building GLINT v${VERSION_CLEAN}"
echo ""

# Step 1: Build GLINT Installer.app
echo "Step 1: Building GLINT Installer.app..."
./build_app.sh
if [ $? -ne 0 ]; then
    echo "❌ Failed to build installer app"
    exit 1
fi
echo ""

# Step 2: Build standalone GLINT.app
echo "Step 2: Building standalone GLINT.app..."
chmod +x build_glint_app.sh
./build_glint_app.sh
if [ $? -ne 0 ]; then
    echo "❌ Failed to build GLINT app"
    exit 1
fi
echo ""

# Step 3: Create installer DMG
echo "Step 3: Creating installer DMG..."
./create_dmg.sh
if [ $? -ne 0 ]; then
    echo "❌ Failed to create installer DMG"
    exit 1
fi
echo ""

# Step 4: Create drag-and-drop DMG
echo "Step 4: Creating drag-and-drop DMG..."
chmod +x create_drag_drop_dmg.sh
./create_drag_drop_dmg.sh
if [ $? -ne 0 ]; then
    echo "❌ Failed to create drag-and-drop DMG"
    exit 1
fi
echo ""

# Move DMGs to parent directory
cd ..

# Rename DMGs to standard names
INSTALLER_DMG=$(ls -t GLINT_Installer_v*.dmg 2>/dev/null | head -1)
DRAGDROP_DMG=$(ls -t GLINT_v*.dmg 2>/dev/null | head -1)

if [ -n "$INSTALLER_DMG" ]; then
    STANDARD_INSTALLER="GLINT_Installer_v${VERSION_CLEAN}.dmg"
    if [ "$INSTALLER_DMG" != "$STANDARD_INSTALLER" ]; then
        mv "$INSTALLER_DMG" "$STANDARD_INSTALLER"
    fi
    echo "✅ Installer DMG: $STANDARD_INSTALLER ($(du -h "$STANDARD_INSTALLER" | cut -f1))"
fi

if [ -n "$DRAGDROP_DMG" ]; then
    STANDARD_DRAGDROP="GLINT_v${VERSION_CLEAN}.dmg"
    if [ "$DRAGDROP_DMG" != "$STANDARD_DRAGDROP" ]; then
        mv "$DRAGDROP_DMG" "$STANDARD_DRAGDROP"
    fi
    echo "✅ Drag-and-drop DMG: $STANDARD_DRAGDROP ($(du -h "$STANDARD_DRAGDROP" | cut -f1))"
fi

echo ""
echo "=========================================="
echo "✅ Build Complete!"
echo "=========================================="
echo ""
echo "Two distribution options created:"
echo ""
echo "1. GLINT_Installer_v${VERSION_CLEAN}.dmg"
echo "   - Full installer with setup wizard"
echo "   - Installs conda environment and dependencies"
echo "   - Recommended for first-time users"
echo ""
echo "2. GLINT_v${VERSION_CLEAN}.dmg"
echo "   - Drag-and-drop installation"
echo "   - Drag GLINT.app to Applications"
echo "   - Requires dependencies already installed"
echo "   - For users who already ran the installer"
echo ""
echo "Distribution strategy:"
echo "  - Provide installer DMG for new users"
echo "  - Provide drag-and-drop DMG for updates"

