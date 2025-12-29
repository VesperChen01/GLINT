#!/bin/bash
# GlueTK Version Update Script
# Usage: ./update_version.sh <new_version>
# Example: ./update_version.sh v0.1.5-beta

if [ -z "$1" ]; then
    echo "Usage: $0 <new_version>"
    echo "Example: $0 v0.1.5-beta"
    echo ""
    echo "Current version:"
    grep -E "__version__" gluetk/_version.py
    exit 1
fi

NEW_VERSION="$1"
OLD_VERSION=$(grep '__version__' gluetk/_version.py | sed 's/.*"\(.*\)".*/\1/')

echo "Updating version: $OLD_VERSION -> $NEW_VERSION"
echo ""

# Update _version.py (main version definition)
sed -i '' "s/__version__ = \".*\"/__version__ = \"$NEW_VERSION\"/" gluetk/_version.py

echo "✅ Updated gluetk/_version.py"

# No longer need to update __init__.py directly, as it will import from _version.py
# No longer need to copy __init__.py to installer bundle or PyMOL startup, as they will import from _version.py

echo ""
echo "Done! New version: $NEW_VERSION"
echo ""
echo "📦 To create new DMG installer, run:"
echo "   ./package_dmg.sh"
echo ""
echo "   This will create: GlueTK_Installer_${NEW_VERSION}.dmg"