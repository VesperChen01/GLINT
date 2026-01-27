#!/bin/bash
# Verification script for GLINT macOS Installer

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR/.."

echo "=========================================="
echo "GLINT macOS Installer - Verification"
echo "=========================================="
echo ""

# Check if DMG exists
DMG_PATH="GLINT_Installer_v0.1.28-beta.dmg"
if [ ! -f "$DMG_PATH" ]; then
    echo "❌ DMG not found: $DMG_PATH"
    exit 1
fi

echo "✅ DMG exists: $DMG_PATH"
echo "   Size: $(du -h "$DMG_PATH" | cut -f1)"
echo ""

# Check if .app exists
APP_PATH="GLINT Installer.app"
if [ ! -d "$APP_PATH" ]; then
    echo "❌ App bundle not found: $APP_PATH"
    exit 1
fi

echo "✅ App bundle exists: $APP_PATH"
echo ""

# Check critical files
echo "Checking critical files..."

FILES=(
    "$APP_PATH/Contents/MacOS/launcher"
    "$APP_PATH/Contents/Resources/GLINT_Installer.py"
    "$APP_PATH/Contents/Resources/glint/__init__.py"
    "$APP_PATH/Contents/Info.plist"
)

ALL_OK=true
for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✅ $file"
    else
        echo "  ❌ Missing: $file"
        ALL_OK=false
    fi
done

echo ""

# Check if launcher is executable
if [ -x "$APP_PATH/Contents/MacOS/launcher" ]; then
    echo "✅ Launcher is executable"
else
    echo "❌ Launcher is not executable"
    ALL_OK=false
fi

echo ""

# Mount DMG and verify contents
echo "Mounting DMG to verify contents..."
MOUNT_OUTPUT=$(hdiutil attach "$DMG_PATH" -readonly -nobrowse 2>&1)
# Extract the mount point - it's after the last tab character on the line with /Volumes
MOUNT_POINT=$(echo "$MOUNT_OUTPUT" | grep "/Volumes" | tail -1 | sed 's/.*[[:space:]]\(\/Volumes\/.*\)$/\1/')

if [ -z "$MOUNT_POINT" ] || [ ! -d "$MOUNT_POINT" ]; then
    echo "❌ Failed to mount DMG or find mount point"
    exit 1
fi

echo "✅ DMG mounted at: $MOUNT_POINT"

# Check if app exists in DMG
if [ -d "$MOUNT_POINT/GLINT Installer.app" ]; then
    echo "✅ App bundle found in DMG"

    # Check if installer script exists in DMG
    if [ -f "$MOUNT_POINT/GLINT Installer.app/Contents/Resources/GLINT_Installer.py" ]; then
        echo "✅ Installer script found in DMG"
    else
        echo "❌ Installer script missing in DMG"
        ALL_OK=false
    fi
else
    echo "❌ App bundle not found in DMG"
    ALL_OK=false
fi

# Unmount DMG
hdiutil detach "$MOUNT_POINT" -quiet
echo "✅ DMG unmounted"

echo ""
echo "=========================================="
if [ "$ALL_OK" = true ]; then
    echo "✅ All checks passed!"
    echo "=========================================="
    echo ""
    echo "The installer is ready for distribution."
    echo ""
    echo "Users can:"
    echo "  1. Download GLINT_Installer_v0.1.28-beta.dmg"
    echo "  2. Open the DMG"
    echo "  3. Double-click 'GLINT Installer.app'"
    echo "  4. Follow the installation wizard"
    exit 0
else
    echo "❌ Some checks failed!"
    echo "=========================================="
    exit 1
fi

