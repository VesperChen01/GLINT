#!/bin/bash
# Build GLINT.app - A standalone app that can be dragged to Applications
# This app will launch PyMOL with GLINT loaded

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR/.."

echo "Building GLINT.app for drag-and-drop installation..."

# Clean old build
rm -rf "GLINT.app"

# Create .app directory structure
APP_NAME="GLINT.app"
CONTENTS_DIR="${APP_NAME}/Contents"
MACOS_DIR="${CONTENTS_DIR}/MacOS"
RESOURCES_DIR="${CONTENTS_DIR}/Resources"

mkdir -p "$MACOS_DIR"
mkdir -p "$RESOURCES_DIR"

echo "Created app bundle structure at: $APP_NAME"

# Copy GLINT source code to Resources
echo "Copying GLINT source code..."
cp -r "glint" "$RESOURCES_DIR/"

# Copy icon
if [ -f "glint/assets/logo.png" ]; then
    cp "glint/assets/logo.png" "$RESOURCES_DIR/AppIcon.png"
    echo "Copied icon"
fi

# Create launcher script that will:
# 1. Find conda
# 2. Activate glint environment
# 3. Launch PyMOL with GLINT
cat > "$MACOS_DIR/GLINT" << 'EOF'
#!/bin/bash
# GLINT Launcher

# Find conda
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
    osascript -e 'display alert "Conda Not Found" message "Please install GLINT first using the GLINT Installer.\n\nOr install Miniconda and run the installer."'
    exit 1
fi

# Check if glint environment exists
eval "$($CONDA_EXE shell.bash hook)"
if ! conda env list | grep -q "^glint "; then
    osascript -e 'display alert "GLINT Not Installed" message "Please run the GLINT Installer first to set up the environment and dependencies."'
    exit 1
fi

# Activate glint environment
conda activate glint

# Set environment variables
export KMP_DUPLICATE_LIB_OK=TRUE
export OMP_NUM_THREADS=1

# Get the Resources directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RESOURCES_DIR="$(dirname "$DIR")/Resources"

# Add GLINT to Python path
export PYTHONPATH="${RESOURCES_DIR}:${PYTHONPATH}"

# Launch PyMOL with GLINT
pymol -d "import sys, os; sys.path.insert(0, '${RESOURCES_DIR}'); import glint; glint.glint_gui()"
EOF

chmod +x "$MACOS_DIR/GLINT"
echo "Created launcher script"

# Create Info.plist
cat > "$CONTENTS_DIR/Info.plist" << 'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleExecutable</key>
    <string>GLINT</string>
    <key>CFBundleIconFile</key>
    <string>AppIcon</string>
    <key>CFBundleIdentifier</key>
    <string>com.glint.app</string>
    <key>CFBundleName</key>
    <string>GLINT</string>
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

# Create icns icon (if sips is available)
if [ -f "$RESOURCES_DIR/AppIcon.png" ]; then
    ICONSET_DIR="$RESOURCES_DIR/AppIcon.iconset"
    mkdir -p "$ICONSET_DIR"
    
    # Generate different sizes
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
    
    # Convert to icns
    iconutil -c icns "$ICONSET_DIR" -o "$RESOURCES_DIR/AppIcon.icns" 2>/dev/null
    rm -rf "$ICONSET_DIR"
    
    if [ -f "$RESOURCES_DIR/AppIcon.icns" ]; then
        echo "Created AppIcon.icns"
    fi
fi

echo ""
echo "✅ GLINT.app created successfully!"
echo "Location: $(pwd)/$APP_NAME"
echo ""
echo "This app can be:"
echo "  1. Dragged to Applications folder"
echo "  2. Launched directly (requires GLINT to be installed first)"
echo ""
echo "Note: Users must run 'GLINT Installer.app' first to install dependencies."

