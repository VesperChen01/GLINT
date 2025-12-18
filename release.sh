#!/bin/bash
# release.sh - GlueTK Release Script
# Combines version update, rebuild installer, and DMG packaging
#
# Usage: ./release.sh <new_version>
# Example: ./release.sh v0.1.6-beta

set -e

if [ -z "$1" ]; then
    CURRENT_VERSION=$(grep '__version__' gluetk/__init__.py | sed 's/.*"\(.*\)".*/\1/')
    echo "GlueTK Release Script"
    echo "====================="
    echo ""
    echo "Current version: $CURRENT_VERSION"
    echo ""
    echo "Usage: $0 <new_version>"
    echo "Example: $0 v0.1.6-beta"
    echo ""
    echo "This script will:"
    echo "  1. Update version in gluetk/__init__.py"
    echo "  2. Update version in Info.plist"
    echo "  3. Rebuild GlueTK Installer.app"
    echo "  4. Create DMG: GlueTK_Installer_<version>.dmg"
    echo "  5. Copy to PyMOL startup directory"
    exit 0
fi

NEW_VERSION="$1"
OLD_VERSION=$(grep '__version__' gluetk/__init__.py | sed 's/.*"\(.*\)".*/\1/')

echo "🚀 GlueTK Release: $OLD_VERSION → $NEW_VERSION"
echo ""

# Step 1: Update version in __init__.py
echo "📝 [1/5] Updating version in source files..."
sed -i '' "s/__version__ = \".*\"/__version__ = \"$NEW_VERSION\"/" gluetk/__init__.py
sed -i '' "s/Version: v[0-9]*\.[0-9]*\.[0-9]*-beta/Version: $NEW_VERSION/" gluetk/__init__.py
sed -i '' "s/GlueTK - Molecular Glue Analyzer v[0-9]*\.[0-9]*\.[0-9]*-beta/GlueTK - Molecular Glue Analyzer $NEW_VERSION/g" gluetk/__init__.py
sed -i '' "s/GlueTK v[0-9]*\.[0-9]*\.[0-9]*-beta/GlueTK $NEW_VERSION/g" gluetk/__init__.py
echo "   ✅ Updated gluetk/__init__.py"

# Step 2: Update version in rebuild_installer.sh (Info.plist template)
echo "📝 [2/5] Updating version in rebuild_installer.sh..."
# Extract version number without 'v' prefix for CFBundleShortVersionString
VERSION_NUM=$(echo "$NEW_VERSION" | sed 's/^v//')
sed -i '' "s/<string>[0-9]*\.[0-9]*\.[0-9]*-beta<\/string>/<string>$VERSION_NUM<\/string>/" rebuild_installer.sh
echo "   ✅ Updated rebuild_installer.sh"

# Step 3: Rebuild Installer.app
echo "🏗️  [3/5] Rebuilding GlueTK Installer.app..."
./rebuild_installer.sh
echo "   ✅ Rebuilt Installer.app"

# Step 4: Package DMG
echo "📦 [4/5] Creating DMG..."
./package_dmg.sh
echo "   ✅ Created GlueTK_Installer_${NEW_VERSION}.dmg"

# Step 5: Copy to PyMOL startup
echo "📋 [5/5] Updating PyMOL startup directory..."
if [ -d ~/.pymol/startup/gluetk ]; then
    rsync -av --exclude='__pycache__' --exclude='*.pyc' --exclude='.DS_Store' \
          gluetk/ ~/.pymol/startup/gluetk/
    echo "   ✅ Updated ~/.pymol/startup/gluetk/"
else
    echo "   ⚠️  PyMOL startup directory not found (will be created on first install)"
fi

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "✅ Release $NEW_VERSION complete!"
echo ""
echo "📦 DMG: GlueTK_Installer_${NEW_VERSION}.dmg"
echo ""
echo "Next steps:"
echo "  1. Test the DMG on a clean system"
echo "  2. Upload to GitHub releases"
echo "  3. Update README if needed"
echo "═══════════════════════════════════════════════════════════"