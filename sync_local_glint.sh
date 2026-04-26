#!/bin/bash
# Sync the working-tree GLINT package into the local PyMOL startup plugin copy.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$SCRIPT_DIR/glint"
TARGET_DIR="${GLINT_LOCAL_TARGET:-$HOME/.pymol/startup/glint}"

if [[ ! -d "$SOURCE_DIR" ]]; then
    echo "Source directory not found: $SOURCE_DIR" >&2
    exit 1
fi

mkdir -p "$(dirname "$TARGET_DIR")"
mkdir -p "$TARGET_DIR"

rsync -av --delete \
    --exclude='__pycache__' \
    --exclude='*.pyc' \
    --exclude='.DS_Store' \
    "$SOURCE_DIR/" "$TARGET_DIR/"

echo "Synced GLINT to $TARGET_DIR"
