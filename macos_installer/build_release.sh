#!/bin/bash
# Complete build and release script for GLINT macOS Installer
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "GLINT macOS Installer - Complete Build"
echo "=========================================="
echo ""

# Step 1: Build the app
echo "Step 1: Building .app bundle..."
./build_app.sh
if [ $? -ne 0 ]; then
    echo "❌ Failed to build app bundle"
    exit 1
fi
echo ""

# Step 2: Create DMG
echo "Step 2: Creating DMG..."
./create_dmg.sh
if [ $? -ne 0 ]; then
    echo "❌ Failed to create DMG"
    exit 1
fi
echo ""

# Step 3: Move DMG to parent directory and rename
cd ..
NEW_DMG=$(ls -t GLINT_Installer_v*.dmg 2>/dev/null | head -1)

if [ -n "$NEW_DMG" ]; then
    # Keep a stable output name for distribution
    STANDARD_NAME="GLINT_Installer_v0.1.28-beta.dmg"
    if [ "$NEW_DMG" != "$STANDARD_NAME" ]; then
        if [ -f "$STANDARD_NAME" ]; then
            echo "Removing old DMG: $STANDARD_NAME"
            rm -f "$STANDARD_NAME"
        fi
        mv "$NEW_DMG" "$STANDARD_NAME"
        echo "Renamed DMG to: $STANDARD_NAME"
    fi
    
    echo ""
    echo "=========================================="
    echo "✅ Build Complete!"
    echo "=========================================="
    echo "DMG Location: $(pwd)/$STANDARD_NAME"
    echo "Size: $(du -h "$STANDARD_NAME" | cut -f1)"
    echo ""
    echo "You can now distribute this DMG to users."
else
    echo "❌ DMG file not found"
    exit 1
fi

