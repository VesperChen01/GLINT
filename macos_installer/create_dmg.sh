#!/bin/bash
# Create DMG installer for GLINT

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR/.."

# Get version from _version.py
VERSION=$(grep "__version__" glint/_version.py | cut -d'"' -f2)
if [ -z "$VERSION" ]; then
    VERSION="0.1.28-beta"
fi

APP_NAME="GLINT Installer.app"
DMG_NAME="GLINT_Installer_v${VERSION}.dmg"
VOLUME_NAME="GLINT Installer"

echo "Creating DMG for GLINT Installer v${VERSION}..."

# Check if app exists
if [ ! -d "$APP_NAME" ]; then
    echo "Error: $APP_NAME not found. Run build_app.sh first."
    exit 1
fi

# Remove old DMG
rm -f "$DMG_NAME"

# Create temporary directory
TMP_DIR=$(mktemp -d)
echo "Using temporary directory: $TMP_DIR"

# Copy app to temp directory
cp -R "$APP_NAME" "$TMP_DIR/"

# Create DMG
echo "Creating DMG..."
hdiutil create -volname "$VOLUME_NAME" \
    -srcfolder "$TMP_DIR" \
    -ov -format UDZO \
    "$DMG_NAME"

# Clean up
rm -rf "$TMP_DIR"

if [ -f "$DMG_NAME" ]; then
    echo ""
    echo "✅ DMG created successfully!"
    echo "Location: $(pwd)/$DMG_NAME"
    echo "Size: $(du -h "$DMG_NAME" | cut -f1)"
else
    echo "❌ Failed to create DMG"
    exit 1
fi

