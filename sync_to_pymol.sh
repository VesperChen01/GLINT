#!/bin/bash
# Sync GLINT changes into the PyMOL startup directory

set -e

PYMOL_STARTUP="$HOME/.pymol/startup/glint"
SOURCE_DIR="./glint"

echo "🔄 Syncing GLINT files into the PyMOL startup directory..."
echo "   Source: $SOURCE_DIR"
echo "   Target: $PYMOL_STARTUP"
echo ""

if [ ! -d "$PYMOL_STARTUP" ]; then
    echo "❌ PyMOL startup directory does not exist: $PYMOL_STARTUP"
    echo "   Please run install_glint.sh first to install GLINT"
    exit 1
fi

if [ ! -d "$SOURCE_DIR" ]; then
    echo "❌ Source directory does not exist: $SOURCE_DIR"
    exit 1
fi

echo "📦 Copying files..."
rsync -av --delete \
      --exclude='__pycache__' \
      --exclude='*.pyc' \
      --exclude='.DS_Store' \
      --exclude='*.sh' \
      "$SOURCE_DIR/" "$PYMOL_STARTUP/"

echo ""
echo "✅ Sync completed."
echo ""
echo "💡 Tips:"
echo "   1. If PyMOL is running, restart it to load updates"
echo "   2. Or run in PyMOL: reinitialize"
echo "   3. Then re-import: import glint; glint.glint_gui()"
echo ""

