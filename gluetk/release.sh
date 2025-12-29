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

# Step 0.5: Verify EC dependencies
echo -e "${BLUE}[0.5/7] Verifying EC dependencies...${NC}"
EC_DEPS_OK=true

if ! command -v pdb2pqr &> /dev/null; then
    echo -e "${YELLOW}⚠️  pdb2pqr not found in PATH${NC}"
    EC_DEPS_OK=false
else
    echo -e "${GREEN}✅ pdb2pqr found: $(command -v pdb2pqr)${NC}"
fi

if ! command -v apbs &> /dev/null; then
    echo -e "${YELLOW}⚠️  apbs not found in PATH${NC}"
    EC_DEPS_OK=false
else
    echo -e "${GREEN}✅ apbs found: $(command -v apbs)${NC}"
fi

if [ "$EC_DEPS_OK" = false ]; then
    echo -e "${YELLOW}⚠️  Some EC dependencies missing. Install with:${NC}"
    echo "   conda install -c conda-forge apbs pdb2pqr"
    echo -e "${YELLOW}Continuing with release (EC analysis will use fallback mode)...${NC}"
fi
echo ""

# Step 1: Update version in __init__.py
echo -e "${BLUE}[1/7] Updating version in source files...${NC}"
sed -i '' "s/__version__ = \".*\"/__version__ = \"$NEW_VERSION\"/" gluetk/__init__.py
sed -i '' "s/Version: v[0-9]*\.[0-9]*\.[0-9]*-beta/Version: $NEW_VERSION/" gluetk/__init__.py
sed -i '' "s/GlueTK - Molecular Glue Analyzer v[0-9]*\.[0-9]*\.[0-9]*-beta/GlueTK - Molecular Glue Analyzer $NEW_VERSION/g" gluetk/__init__.py
sed -i '' "s/GlueTK v[0-9]*\.[0-9]*\.[0-9]*-beta/GlueTK $NEW_VERSION/g" gluetk/__init__.py
echo "   ✅ Updated gluetk/__init__.py"

# Step 2: Update version in rebuild_installer.sh (Info.plist template)
echo -e "${BLUE}[2/7] Updating version in rebuild_installer.sh...${NC}"
# Extract version number without 'v' prefix for CFBundleShortVersionString
VERSION_NUM=$(echo "$NEW_VERSION" | sed 's/^v//')
sed -i '' "s/<string>[0-9]*\.[0-9]*\.[0-9]*-beta<\/string>/<string>$VERSION_NUM<\/string>/" rebuild_installer.sh
echo -e "${GREEN}   ✅ Updated rebuild_installer.sh${NC}"

# Step 3: Rebuild Installer.app
echo -e "${BLUE}[3/7] Rebuilding GlueTK Installer.app...${NC}"
./rebuild_installer.sh
echo -e "${GREEN}   ✅ Rebuilt Installer.app${NC}"

# Step 4: Package DMG
echo -e "${BLUE}[4/7] Creating DMG...${NC}"
./package_dmg.sh
echo -e "${GREEN}   ✅ Created GlueTK_Installer_${NEW_VERSION}.dmg${NC}"

# Step 5: Copy to PyMOL startup
echo -e "${BLUE}[5/7] Updating PyMOL startup directory...${NC}"
if [ -d ~/.pymol/startup/gluetk ]; then
    rsync -av --exclude='__pycache__' --exclude='*.pyc' --exclude='.DS_Store' \
          gluetk/ ~/.pymol/startup/gluetk/
    echo -e "${GREEN}   ✅ Updated ~/.pymol/startup/gluetk/${NC}"
else
    echo -e "${YELLOW}   ⚠️  PyMOL startup directory not found (will be created on first install)${NC}"
fi

# Step 6: Install EC dependencies (optional)
echo -e "${BLUE}[6/7] Installing EC dependencies...${NC}"
if [ -f "install_ec_dependencies.sh" ]; then
    bash install_ec_dependencies.sh
    echo -e "${GREEN}   ✅ EC dependencies installed${NC}"
else
    echo -e "${YELLOW}   ⚠️  install_ec_dependencies.sh not found${NC}"
fi

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