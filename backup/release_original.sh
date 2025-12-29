#!/bin/bash
# release.sh - GlueTK Release Script
# Combines version update, rebuild installer, and DMG packaging
#
# Usage: ./release.sh <new_version>
# Example: ./release.sh v0.1.6-beta

set -e

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check conda environment
echo -e "${BLUE}[0/7] Checking conda environment...${NC}"
if [ -z "$CONDA_PREFIX" ]; then
    echo -e "${YELLOW}⚠️  Conda environment not activated${NC}"
    echo "Please activate the gluetk environment:"
    echo "  conda activate gluetk"
    exit 1
fi
echo -e "${GREEN}✅ Conda environment: $CONDA_PREFIX${NC}"

if [ -z "$1" ]; then
    CURRENT_VERSION=$(python3 -c "from gluetk._version import __version__; print(__version__)")
    echo "GlueTK Release Script"
    echo "====================="
    echo ""
    echo "Current version: $CURRENT_VERSION"
    echo ""
    echo "Usage: $0 <new_version>"
    echo "Example: $0 v0.1.6-beta"
    echo ""
    echo "This script will:"
    echo "  1. Update version in gluetk/_version.py"
    echo "  2. Rebuild GlueTK Installer.app"
    echo "  3. Create DMG: GlueTK_Installer_<version>.dmg"
    echo "  4. Copy to PyMOL startup directory"
    exit 0
fi

NEW_VERSION="$1"
OLD_VERSION=$(python3 -c "from gluetk._version import __version__; print(__version__)")

echo "🚀 GlueTK Release: $OLD_VERSION → $NEW_VERSION"
echo ""

# Step 1: Update version in _version.py
echo -e "${BLUE}[1/5] Updating version in source files...${NC}"
# Update _version.py (main version definition)
sed -i '' "s/__version__ = \".*\"/__version__ = \"$NEW_VERSION\"/" gluetk/_version.py
echo "   ✅ Updated gluetk/_version.py"

# Step 2: Rebuild Installer.app
echo -e "${BLUE}[2/5] Rebuilding GlueTK Installer.app...${NC}"

# Step 4: Copy to PyMOL startup
echo -e "${BLUE}[4/5] Updating PyMOL startup directory...${NC}"
if [ -d ~/.pymol/startup/gluetk ]; then
    rsync -av --exclude='__pycache__' --exclude='*.pyc' --exclude='.DS_Store' \
          gluetk/ ~/.pymol/startup/gluetk/
    echo -e "${GREEN}   ✅ Updated ~/.pymol/startup/gluetk/${NC}"
else
    echo -e "${YELLOW}   ⚠️  PyMOL startup directory not found (will be created on first install)${NC}"
fi

# Step 5: Clean up (optional)
echo -e "${BLUE}[5/5] Cleaning up...${NC}"
rm -rf build_dmg/GlueTK\ Installer.app
echo -e "${GREEN}   ✅ Cleaned up build artifacts${NC}"

echo ""
echo -e "${BLUE}═══════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}✅ Release $NEW_VERSION complete!${NC}"
echo ""
echo -e "${BLUE}📦 DMG: GlueTK_Installer_${NEW_VERSION}.dmg${NC}"
echo ""
echo -e "${BLUE}Next steps:${NC}"
echo "  1. Test the DMG on a clean system"
echo "  2. Upload to GitHub releases"
echo "  3. Update README if needed"
echo -e "${BLUE}═══════════════════════════════════════════════════════════${NC}"