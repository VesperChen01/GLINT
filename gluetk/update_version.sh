#!/bin/bash
# GlueTK Version Update Script
# Usage: ./update_version.sh <new_version>
# Example: ./update_version.sh v0.1.5-beta

if [ -z "$1" ]; then
    echo "Usage: $0 <new_version>"
    echo "Example: $0 v0.1.5-beta"
    echo ""
    echo "Current version:"
    grep -E "__version__|v0\.[0-9]+\.[0-9]+-beta" gluetk/__init__.py | head -2
    exit 1
fi

NEW_VERSION="$1"
OLD_VERSION=$(grep '__version__' gluetk/__init__.py | sed 's/.*"\(.*\)".*/\1/')

echo "Updating version: $OLD_VERSION -> $NEW_VERSION"
echo ""

# Update __init__.py (main version definition)
sed -i '' "s/__version__ = \".*\"/__version__ = \"$NEW_VERSION\"/" gluetk/__init__.py
sed -i '' "s/Version: v[0-9]*\.[0-9]*\.[0-9]*-beta/Version: $NEW_VERSION/" gluetk/__init__.py

# Update banner messages in __init__.py
sed -i '' "s/GlueTK - Molecular Glue Analyzer v[0-9]*\.[0-9]*\.[0-9]*-beta/GlueTK - Molecular Glue Analyzer $NEW_VERSION/g" gluetk/__init__.py
sed -i '' "s/GlueTK v[0-9]*\.[0-9]*\.[0-9]*-beta/GlueTK $NEW_VERSION/g" gluetk/__init__.py

echo "✅ Updated gluetk/__init__.py"

# Copy to installer and PyMOL startup
cp gluetk/__init__.py "build_dmg/GlueTK Installer.app/Contents/Resources/gluetk/__init__.py" 2>/dev/null && echo "✅ Updated installer bundle"
cp gluetk/__init__.py ~/.pymol/startup/gluetk/__init__.py 2>/dev/null && echo "✅ Updated PyMOL startup"

echo ""
echo "Done! New version: $NEW_VERSION"
echo ""
echo "📦 To create new DMG installer, run:"
echo "   ./package_dmg.sh"
echo ""
echo "   This will create: GlueTK_Installer_${NEW_VERSION}.dmg"