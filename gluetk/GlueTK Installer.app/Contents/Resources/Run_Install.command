#!/bin/bash
# Optional: Terminal-based GlueTK installer entry point (fallback)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GLUETK_DIR="${SCRIPT_DIR}/gluetk"

if [ ! -d "$GLUETK_DIR" ]; then
  echo "[GlueTK Installer] ERROR: GlueTK sources not found in: $GLUETK_DIR"
  echo "If you are running from a source checkout, please run: bash gluetk/install.sh"
  read -n1 -p "Press any key to exit..." _
  exit 1
fi

cd "$GLUETK_DIR"

if [ -x "install.sh" ]; then
  bash install.sh
else
  echo "[GlueTK Installer] ERROR: install.sh not found or not executable."
  read -n1 -p "Press any key to exit..." _
  exit 1
fi

read -n1 -p "Press any key to close this window..." _
