#!/bin/bash
# Create a drag-and-drop DMG installer for GLINT
# This creates a DMG with the app and a link to /Applications
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR/.."

# Get version from source tree
VERSION=""
if [ -f "glint/_version.py" ]; then
    VERSION=$(grep "__version__" glint/_version.py | cut -d'"' -f2 || true)
fi
if [ -z "$VERSION" ]; then
    VERSION="0.2.3"
fi
VERSION_CLEAN="${VERSION#v}"

APP_NAME="GLINT.app"
DMG_NAME="GLINT_v${VERSION_CLEAN}.dmg"
VOLUME_NAME="GLINT"
SOURCE_APP="GLINT.app"

echo "Creating drag-and-drop DMG for GLINT v${VERSION_CLEAN}..."

# Check if source app exists
if [ ! -d "$SOURCE_APP" ]; then
    echo "Error: $SOURCE_APP not found. Run build_glint_app.sh first."
    exit 1
fi

# Remove old DMG
rm -f "$DMG_NAME"

# Create temporary directory for DMG contents
TMP_DIR=$(mktemp -d)
echo "Using temporary directory: $TMP_DIR"

# Copy app to temp directory and rename to GLINT.app
cp -R "$SOURCE_APP" "$TMP_DIR/$APP_NAME"

# Create symbolic link to /Applications
ln -s /Applications "$TMP_DIR/Applications"

# Create a .background folder for custom background (optional)
mkdir -p "$TMP_DIR/.background"

# If you have a custom background image, copy it here
# cp path/to/background.png "$TMP_DIR/.background/background.png"

# Create temporary DMG
TEMP_DMG="temp_${DMG_NAME}"
echo "Creating temporary DMG..."
hdiutil create -volname "$VOLUME_NAME" \
    -srcfolder "$TMP_DIR" \
    -ov -format UDRW \
    "$TEMP_DMG"

# Mount the temporary DMG (use fixed mountpoint to avoid parsing issues)
echo "Mounting temporary DMG..."
MOUNT_DIR="$TMP_DIR/mount"
mkdir -p "$MOUNT_DIR"
hdiutil attach "$TEMP_DMG" -mountpoint "$MOUNT_DIR" -nobrowse >/dev/null

echo "Mounted at: $MOUNT_DIR"

# Set custom view options using AppleScript (best effort)
echo "Setting DMG window properties..."
if ! osascript <<EOF
with timeout of 60 seconds
    tell application "Finder"
        tell disk "$VOLUME_NAME"
            open
            set current view of container window to icon view
            set toolbar visible of container window to false
            set statusbar visible of container window to false
            set the bounds of container window to {100, 100, 700, 500}
            set viewOptions to the icon view options of container window
            set arrangement of viewOptions to not arranged
            set icon size of viewOptions to 128

            -- Set background only if the file exists
            try
                set background picture of viewOptions to file ".background:background.png"
            end try

            -- Position the app icon
            set position of item "$APP_NAME" of container window to {150, 200}

            -- Position the Applications link
            set position of item "Applications" of container window to {450, 200}

            close
            open
            update without registering applications
            delay 2
        end tell
    end tell
end timeout
EOF
then
    echo "⚠️ Finder customization timed out; continuing with default DMG layout."
fi

# Sync and unmount
sync
echo "Unmounting temporary DMG..."
hdiutil detach "$MOUNT_DIR" -quiet

# Convert to compressed read-only DMG
echo "Converting to final DMG..."
hdiutil convert "$TEMP_DMG" -format UDZO -o "$DMG_NAME"

# Clean up
rm -rf "$TMP_DIR"
rm -f "$TEMP_DMG"

if [ -f "$DMG_NAME" ]; then
    echo ""
    echo "✅ Drag-and-drop DMG created successfully!"
    echo "Location: $(pwd)/$DMG_NAME"
    echo "Size: $(du -h "$DMG_NAME" | cut -f1)"
    echo ""
    echo "Users can now:"
    echo "  1. Open the DMG"
    echo "  2. Drag GLINT.app to Applications folder"
    echo "  3. Launch GLINT from Applications"
else
    echo "❌ Failed to create DMG"
    exit 1
fi

