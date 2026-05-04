#!/bin/bash
# GLINT one-click installer
# Supports multiple PyMOL installation layouts with more robust dependency checks

set -e

# Configuration
ENV_NAME="glint"
PYTHON_VERSION="3.11"
INSTALL_DIR="$HOME/.pymol/startup/glint"
DESKTOP_APP="$HOME/Desktop/GLINT.app"

# Terminal colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🧬 GLINT installer${NC}"
echo "================================"

# 0. Detect operating system and architecture
echo -e "\n${BLUE}[0/6]${NC} Detecting system environment and architecture..."
OS_TYPE=$(uname -s)
ARCH_TYPE=$(uname -m)

echo "   Operating system: $OS_TYPE"
echo "   Architecture: $ARCH_TYPE"

if [ "$OS_TYPE" = "Darwin" ]; then
    if [ "$ARCH_TYPE" = "arm64" ]; then
        echo -e "${GREEN}✅ Apple Silicon (M1/M2/M3) detected${NC}"
        HOMEBREW_PREFIX="/opt/homebrew"
    else
        echo -e "${GREEN}✅ Intel Mac detected${NC}"
        HOMEBREW_PREFIX="/usr/local"
    fi

    # Add Homebrew to PATH so brew can be discovered reliably later
    if [ -d "$HOMEBREW_PREFIX/bin" ]; then
        export PATH="$HOMEBREW_PREFIX/bin:$PATH"
        echo -e "${GREEN}✅ Configured Homebrew path: $HOMEBREW_PREFIX/bin${NC}"
    fi
fi


# 1. Check Conda
echo -e "\n${BLUE}[1/6]${NC} Checking Conda environment..."
CONDA_EXE=""
if command -v conda &> /dev/null; then
    CONDA_EXE=$(command -v conda)
else
    # Try common install paths
    for p in "$HOME/miniconda3/bin/conda" \
             "$HOME/anaconda3/bin/conda" \
             "/opt/miniconda3/bin/conda" \
             "/opt/anaconda3/bin/conda" \
             "/usr/local/bin/conda" \
             "/opt/homebrew/bin/conda" \
             "/opt/homebrew/Caskroom/miniconda/base/bin/conda" \
             "$HOME/opt/miniconda3/bin/conda"; do
        if [ -x "$p" ]; then
            CONDA_EXE="$p"
            break
        fi
    done
fi

if [ -z "$CONDA_EXE" ]; then
    echo -e "${RED}❌ Conda was not found${NC}"
    echo "Please install Miniconda first: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo -e "${GREEN}✅ Found Conda: $CONDA_EXE${NC}"

# Initialize conda
eval "$($CONDA_EXE shell.bash hook)"

# 2. Check or create the conda environment
echo -e "\n${BLUE}[2/6]${NC} Checking Conda environment '$ENV_NAME'..."
if conda env list | grep -q "^${ENV_NAME} "; then
    echo -e "${YELLOW}⚠️  Environment already exists; dependencies will be updated${NC}"
else
    echo "   Creating a new environment..."
    conda create -n "$ENV_NAME" python=$PYTHON_VERSION -y
    echo -e "${GREEN}✅ Environment created successfully${NC}"
fi

# 3. Detect PyMOL installation mode
echo -e "\n${BLUE}[3/6]${NC} Detecting PyMOL..."
PYMOL_TYPE=""
PYMOL_PATH=""
HAS_SYSTEM_PYMOL=false

# 3a. Check system PyMOL.app (macOS commercial/open-source app bundle)
if [ -d "/Applications/PyMOL.app" ]; then
    PYMOL_APP_BIN="/Applications/PyMOL.app/Contents/MacOS/PyMOL"
    if [ -x "$PYMOL_APP_BIN" ]; then
        HAS_SYSTEM_PYMOL=true
        echo -e "${GREEN}✅ Detected system PyMOL.app: $PYMOL_APP_BIN${NC}"
    fi
fi

# 3b. Prefer installing PyMOL into the conda environment to avoid dependency conflicts
echo "   Installing PyMOL into the conda environment (recommended)..."
echo "   This may take a few minutes..."

# Try installing pymol-open-source
if conda install -n "$ENV_NAME" -c conda-forge pymol-open-source -y; then
    # Verify installation
    if conda run -n "$ENV_NAME" which pymol &> /dev/null; then
        PYMOL_TYPE="conda"
        PYMOL_PATH="conda"
        echo -e "${GREEN}✅ PyMOL installed into the conda environment${NC}"

        if [ "$HAS_SYSTEM_PYMOL" = true ]; then
            echo -e "${YELLOW}💡 Note: a system PyMOL.app was also detected, but GLINT will use the conda-managed PyMOL${NC}"
            echo -e "${YELLOW}   This avoids dependency conflicts and keeps all features available${NC}"
        fi
    else
        echo -e "${RED}❌ PyMOL installation verification failed${NC}"

        if [ "$HAS_SYSTEM_PYMOL" = true ]; then
            echo -e "${YELLOW}⚠️  Falling back to system PyMOL.app (dependency conflicts are possible)${NC}"
            PYMOL_TYPE="app"
            PYMOL_PATH="$PYMOL_APP_BIN"
        else
            echo -e "${RED}❌ Unable to install PyMOL${NC}"
            echo "Please install PyMOL.app manually: https://pymol.org/"
            exit 1
        fi
    fi
else
    echo -e "${YELLOW}⚠️  Could not install PyMOL via conda${NC}"

    if [ "$HAS_SYSTEM_PYMOL" = true ]; then
        echo -e "${YELLOW}   System PyMOL.app will be used instead (dependency conflicts are possible)${NC}"
        PYMOL_TYPE="app"
        PYMOL_PATH="$PYMOL_APP_BIN"
    else
        echo -e "${RED}❌ No usable PyMOL installation was found${NC}"
        echo "Please install PyMOL.app manually into /Applications/"
        exit 1
    fi
fi

# 4. Install dependencies
echo -e "\n${BLUE}[4/6]${NC} Installing dependencies..."
echo "   This may take a few minutes. Please wait..."

# Core dependencies installed through conda
CONDA_PACKAGES=(
    "rdkit"
    "scipy"
    "matplotlib"
    "pillow"
    "numpy<2.0.0" # Prevent incompatibilities with legacy packages on numpy v2
    "pandas"
    "seaborn"
    "pyqt"
    "openbabel"
    "requests"
    "vina"
    "pdb2pqr"
    "apbs"  # Required for electrostatic complementarity analysis
    "open3d"  # Preferred surface backend across macOS/Windows/Linux
)

# Temporarily suppress InsecureRequestWarning during installation
export PYTHONWARNINGS="ignore::urllib3.exceptions.InsecureRequestWarning"
conda install -n "$ENV_NAME" -c conda-forge "${CONDA_PACKAGES[@]}" python=$PYTHON_VERSION -y
unset PYTHONWARNINGS # Clear the override after installation

# Surface-analysis extras
echo "   Installing optional surface-analysis extras (scikit-image)..."
conda run -n "$ENV_NAME" python -m pip install scikit-image --quiet --disable-pip-version-check || {
    echo -e "${YELLOW}⚠️  scikit-image installation failed; GLINT will use the built-in EDTSurf fallback path${NC}"
}

echo -e "${GREEN}✅ Dependency installation completed${NC}"

# 4a. Detect macOS architecture and configure Homebrew paths
echo -e "\n${BLUE}[4a]${NC} Detecting macOS architecture and configuring Homebrew paths..."
MACOS_OS_TYPE=$(uname -s)
MACOS_ARCH_TYPE=$(uname -m)
HOMEBREW_PREFIX="/usr/local"
HOMEBREW_BIN=""
HOMEBREW_SBIN=""
BREW_CMD=""

echo "   Detected system: $MACOS_OS_TYPE"
echo "   Detected architecture: $MACOS_ARCH_TYPE"

if [ "$MACOS_OS_TYPE" = "Darwin" ]; then
    if [ "$MACOS_ARCH_TYPE" = "arm64" ]; then
        HOMEBREW_PREFIX="/opt/homebrew"
        echo -e "${GREEN}✅ Apple Silicon macOS detected${NC}"
    elif [ "$MACOS_ARCH_TYPE" = "x86_64" ]; then
        HOMEBREW_PREFIX="/usr/local"
        echo -e "${GREEN}✅ Intel macOS detected${NC}"
    else
        HOMEBREW_PREFIX="/usr/local"
        echo -e "${YELLOW}⚠️  Unrecognized macOS architecture: $MACOS_ARCH_TYPE; defaulting to Homebrew prefix $HOMEBREW_PREFIX${NC}"
    fi
else
    HOMEBREW_PREFIX="/usr/local"
    echo -e "${YELLOW}⚠️  Current system is not macOS; defaulting to Homebrew prefix $HOMEBREW_PREFIX${NC}"
fi

HOMEBREW_BIN="$HOMEBREW_PREFIX/bin"
HOMEBREW_SBIN="$HOMEBREW_PREFIX/sbin"
export HOMEBREW_BIN
export HOMEBREW_SBIN
export PATH="$HOMEBREW_BIN:$HOMEBREW_SBIN:$PATH"

echo "   Homebrew prefix: $HOMEBREW_PREFIX"
echo "   Homebrew bin path: $HOMEBREW_BIN"
echo "   Homebrew sbin path: $HOMEBREW_SBIN"

if [ -x "$HOMEBREW_BIN/brew" ]; then
    BREW_CMD="$HOMEBREW_BIN/brew"
elif command -v brew &> /dev/null; then
    BREW_CMD=$(command -v brew)
else
    BREW_CMD=""
fi

if [ -n "$BREW_CMD" ]; then
    echo -e "${GREEN}✅ brew command found: $BREW_CMD${NC}"
else
    echo -e "${YELLOW}⚠️  brew command not found${NC}"
fi

# 4b. Check and install GCC (required by HADDOCK3/CNS)
echo -e "\n${BLUE}[4b]${NC} Checking GCC dependency..."
if [ -z "$BREW_CMD" ]; then
    echo -e "${YELLOW}⚠️  Homebrew was not found; skipping automatic GCC installation, and HADDOCK3/CNS may not work correctly${NC}"
else
    if ! "$BREW_CMD" list gcc &> /dev/null; then
        echo "   Installing GCC (required by HADDOCK3)..."
        echo "   ${YELLOW}This may take 5-10 minutes. Please wait...${NC}"
        if "$BREW_CMD" install gcc; then
            echo -e "${GREEN}✅ GCC installed successfully${NC}"
        else
            echo -e "${YELLOW}⚠️  GCC installation failed; HADDOCK3 may not work correctly${NC}"
        fi
    else
        echo -e "${GREEN}✅ GCC is already installed${NC}"
    fi
fi

echo "   GCC bin directory for the current architecture: $HOMEBREW_BIN"
GCC_BIN_DIR="$HOMEBREW_BIN"
LATEST_GCC=$(ls "$GCC_BIN_DIR"/gcc-[0-9]* 2>/dev/null | sort -V | tail -n 1)

export CPPFLAGS="-I$HOMEBREW_PREFIX/include ${CPPFLAGS:-}"
export LDFLAGS="-L$HOMEBREW_PREFIX/lib ${LDFLAGS:-}"

if [ -n "$LATEST_GCC" ] && [ -x "$LATEST_GCC" ]; then
    GCC_VERSION_SUFFIX=$(basename "$LATEST_GCC")
    GCC_VERSION_SUFFIX="${GCC_VERSION_SUFFIX#gcc-}"

    export CC="$LATEST_GCC"

    if [ -x "$GCC_BIN_DIR/g++-$GCC_VERSION_SUFFIX" ]; then
        export CXX="$GCC_BIN_DIR/g++-$GCC_VERSION_SUFFIX"
    else
        unset CXX
    fi

    if [ -x "$GCC_BIN_DIR/gfortran-$GCC_VERSION_SUFFIX" ]; then
        export FC="$GCC_BIN_DIR/gfortran-$GCC_VERSION_SUFFIX"
    else
        unset FC
    fi

    echo -e "${GREEN}✅ Found versioned GCC: $LATEST_GCC${NC}"
    echo "   CC=$CC"
    echo "   CXX=${CXX:-unset}"
    echo "   FC=${FC:-unset}"
    echo "   CPPFLAGS=$CPPFLAGS"
    echo "   LDFLAGS=$LDFLAGS"
else
    echo -e "${YELLOW}⚠️  No versioned gcc executable found: $GCC_BIN_DIR/gcc-[0-9]*${NC}"
    echo "   CC=${CC:-unset}"
    echo "   CXX=${CXX:-unset}"
    echo "   FC=${FC:-unset}"
    echo "   CPPFLAGS=$CPPFLAGS"
    echo "   LDFLAGS=$LDFLAGS"
fi

# 4c. Install HADDOCK3 via pip (protein-protein docking engine)
echo -e "\n${BLUE}[4c]${NC} Installing HADDOCK3..."
echo "   ${YELLOW}Installing HADDOCK3 with pip...${NC}"
if conda run -n "$ENV_NAME" python -m pip install -U haddock3 --quiet --disable-pip-version-check; then
    echo -e "${GREEN}✅ HADDOCK3 (pip) installed successfully${NC}"

    # Repair and validate haddock3 permissions
    HADDOCK3_DIR=""
    CONDA_BIN_DIR="$(conda info --base)/envs/$ENV_NAME/bin"
    if [ "$OS_TYPE" = "Darwin" ] && [ -n "$HOMEBREW_PREFIX" ] && ls "$HOMEBREW_PREFIX/bin/haddock3"* &> /dev/null; then
        HADDOCK3_DIR="$HOMEBREW_PREFIX/bin"
    elif [ -d "$CONDA_BIN_DIR" ] && ls "$CONDA_BIN_DIR/haddock3"* &> /dev/null; then
        HADDOCK3_DIR="$CONDA_BIN_DIR"
    fi

    if [ -n "$HADDOCK3_DIR" ]; then
        find "$HADDOCK3_DIR" -maxdepth 1 -name "haddock3*" -type f -exec chmod +x {} \;
        echo -e "${GREEN}✅ Executable permissions updated: $HADDOCK3_DIR/haddock3*${NC}"
    else
        echo -e "${YELLOW}⚠️  No haddock3 executables found; you may need to grant permissions manually${NC}"
    fi

    # Validate the haddock3 command
    if conda run -n "$ENV_NAME" haddock3 -h &> /dev/null; then
        echo -e "${GREEN}✅ haddock3 command validation passed${NC}"
    else
        echo -e "${YELLOW}⚠️  haddock3 command is unavailable; please verify the installation${NC}"
    fi

else
    echo -e "${YELLOW}⚠️  HADDOCK3 installation failed; protein-protein docking features will be unavailable${NC}"
    echo "   You can retry later with: conda run -n $ENV_NAME python -m pip install -U haddock3"
fi

# 4d. Verify surface/EC analysis dependencies
echo -e "\n${BLUE}[4d]${NC} Verifying surface/EC analysis dependencies..."
EC_DEPS_OK=true

# Check Open3D
if conda run -n "$ENV_NAME" python -c "import open3d" 2>/dev/null; then
    echo -e "${GREEN}✅ Open3D is available${NC}"
else
    echo -e "${YELLOW}⚠️  Open3D is unavailable; surface analysis will use built-in EDTSurf${NC}"
fi

# Check PDB2PQR
if conda run -n "$ENV_NAME" python -c "import pdb2pqr" 2>/dev/null; then
    echo -e "${GREEN}✅ PDB2PQR Python API is available${NC}"
else
    echo -e "${YELLOW}⚠️  PDB2PQR Python API is unavailable${NC}"
    EC_DEPS_OK=false
fi

if command -v pdb2pqr &> /dev/null; then
    echo -e "${GREEN}✅ PDB2PQR CLI is available${NC}"
else
    echo -e "${YELLOW}⚠️  PDB2PQR CLI is unavailable${NC}"
fi

# Check APBS
if conda run -n "$ENV_NAME" python -c "from apbs_binary import run_apbs" 2>/dev/null; then
    echo -e "${GREEN}✅ APBS (apbs-binary) is available${NC}"
elif conda run -n "$ENV_NAME" python -c "import apbs" 2>/dev/null; then
    echo -e "${GREEN}✅ APBS (Python API) is available${NC}"
elif command -v apbs &> /dev/null; then
    echo -e "${GREEN}✅ APBS CLI is available${NC}"
else
    echo -e "${YELLOW}⚠️  APBS is unavailable; EC analysis will not run${NC}"
    EC_DEPS_OK=false
fi

if [ "$EC_DEPS_OK" = false ]; then
    echo -e "${YELLOW}💡 Tip: run 'bash install_ec_dependencies.sh' later to finish EC analysis setup${NC}"
fi

# 5. Copy GLINT files
echo -e "\n${BLUE}[5/6]${NC} Installing GLINT plugin..."
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SOURCE_DIR="$SCRIPT_DIR/glint"

if [ ! -d "$SOURCE_DIR" ]; then
    echo -e "${RED}❌ GLINT source directory not found: $SOURCE_DIR${NC}"
    exit 1
fi

if [ ! -f "$SOURCE_DIR/__init__.py" ]; then
    echo -e "${RED}❌ GLINT source files not found${NC}"
    exit 1
fi

# Create installation directory
mkdir -p "$(dirname "$INSTALL_DIR")"

# Copy files
echo "   Copying files to $INSTALL_DIR ..."
if [ -d "$INSTALL_DIR" ]; then
    rm -rf "$INSTALL_DIR"
fi

mkdir -p "$INSTALL_DIR"
rsync -av --exclude='__pycache__' \
          --exclude='*.pyc' \
          --exclude='.DS_Store' \
          --exclude='.git' \
          --exclude='*.sh' \
          "$SOURCE_DIR/" "$INSTALL_DIR/"

echo -e "${GREEN}✅ GLINT files installed${NC}"

# 5b. Create PyMOL startup script
echo -e "\n${BLUE}[5b]${NC} Configuring PyMOL startup script..."
PYMOL_STARTUP_DIR="$HOME/.pymol/startup"
mkdir -p "$PYMOL_STARTUP_DIR"

# Copy startup script
cp "$SCRIPT_DIR/pymol_startup_glint.py" "$PYMOL_STARTUP_DIR/01_glint.py"
echo -e "${GREEN}✅ PyMOL startup script configured${NC}"

# 6. Create launcher (based on PyMOL mode)
echo -e "\n${BLUE}[6/6]${NC} Creating desktop application..."

if [ -d "$DESKTOP_APP" ]; then
    rm -rf "$DESKTOP_APP"
fi

CONTENTS="$DESKTOP_APP/Contents"
MACOS="$CONTENTS/MacOS"
RESOURCES="$CONTENTS/Resources"

mkdir -p "$MACOS"
mkdir -p "$RESOURCES"

# Copy icon
if [ -f "$SOURCE_DIR/assets/AppIcon.icns" ]; then
    cp "$SOURCE_DIR/assets/AppIcon.icns" "$RESOURCES/AppIcon.icns"
fi

# Create Info.plist
cat > "$CONTENTS/Info.plist" << 'PLIST_EOF'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleExecutable</key>
    <string>launcher</string>
    <key>CFBundleIconFile</key>
    <string>AppIcon</string>
    <key>CFBundleIdentifier</key>
    <string>com.vesper.glint</string>
    <key>CFBundleName</key>
    <string>GLINT</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>CFBundleShortVersionString</key>
    <string>0.3.2</string>
</dict>
</plist>
PLIST_EOF


# Create launcher script (based on PyMOL mode)
cat > "$MACOS/launcher" << 'LAUNCHER_EOF'
#!/bin/bash
# GLINT Launcher (macOS .app)

CONDA_EXE=""
if command -v conda &> /dev/null; then
    CONDA_EXE=$(command -v conda)
else
    for p in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" "/opt/miniconda3/bin/conda" "/opt/anaconda3/bin/conda" "/usr/local/bin/conda" "/opt/homebrew/bin/conda" "/opt/homebrew/Caskroom/miniconda/base/bin/conda" "$HOME/opt/miniconda3/bin/conda"; do
        if [ -x "$p" ]; then
            CONDA_EXE="$p"
            break
        fi
    done
fi

if [ -z "$CONDA_EXE" ]; then
    osascript -e 'display alert "GLINT Error" message "Conda not found.\n\nPlease install Miniconda first:\nhttps://docs.conda.io/en/latest/miniconda.html"'
    exit 1
fi

CONDA_BASE=$("$CONDA_EXE" info --base 2>/dev/null)
ENV_PATH="$CONDA_BASE/envs/glint"

if [ ! -d "$ENV_PATH" ]; then
    osascript -e 'display alert "GLINT Error" message "Conda environment glint not found.\n\nPlease run the GLINT installer first."'
    exit 1
fi

export KMP_DUPLICATE_LIB_OK=TRUE
export OMP_NUM_THREADS=1
export QT_MAC_WANTS_LAYER=1

eval "$($CONDA_EXE shell.bash hook)" 2>/dev/null

conda activate glint 2>/dev/null || {
    osascript -e 'display alert "GLINT Error" message "Failed to activate conda environment glint.\n\nPlease check your conda installation."'
    exit 1
}

if [ -n "$CONDA_PREFIX" ]; then
    CONDA_SP=$(python -c "import site; print(site.getsitepackages()[0])" 2>/dev/null)
    if [ -n "$CONDA_SP" ] && [ -d "$CONDA_SP" ]; then
        export PYTHONPATH="${CONDA_SP}:${PYTHONPATH}"
    else
        for pyver in 3.12 3.11 3.10 3.9; do
            SP="${CONDA_PREFIX}/lib/python${pyver}/site-packages"
            if [ -d "$SP" ]; then
                export PYTHONPATH="${SP}:${PYTHONPATH}"
                break
            fi
        done
    fi
fi

LAUNCHER_EOF

# Append launch command based on the detected PyMOL type
if [ "$PYMOL_TYPE" = "app" ]; then
    cat >> "$MACOS/launcher" << 'LAUNCHER_EOF'
/Applications/PyMOL.app/Contents/MacOS/PyMOL -d "
import sys, os
startup_path = os.path.expanduser('~/.pymol/startup')
if startup_path not in sys.path:
    sys.path.insert(0, startup_path)

try:
    import glint
    glint.glint_gui()
except Exception:
    raise
"
LAUNCHER_EOF
else
    cat >> "$MACOS/launcher" << 'LAUNCHER_EOF'
pymol -d "
import sys, os
startup_path = os.path.expanduser('~/.pymol/startup')
if startup_path not in sys.path:
    sys.path.insert(0, startup_path)
try:
    import glint
    glint.glint_gui()
except Exception:
    raise
"
LAUNCHER_EOF
fi

chmod +x "$MACOS/launcher"

echo -e "${GREEN}✅ Desktop application created: $DESKTOP_APP${NC}"

# Done
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}🎉 GLINT installation completed successfully${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "How to launch:"
echo "  1. Double-click GLINT.app on your Desktop"
echo "  2. Or run inside PyMOL:"
echo "     import glint"
echo "     glint.glint_gui()"
echo ""
echo "PyMOL mode: $PYMOL_TYPE"
if [ "$PYMOL_TYPE" = "conda" ]; then
    echo -e "${GREEN}✅ Using the conda-managed PyMOL (recommended)${NC}"
    echo -e "${GREEN}   Dependencies are isolated to avoid conflicts with a system PyMOL${NC}"
    if [ "$HAS_SYSTEM_PYMOL" = true ]; then
        echo -e "${BLUE}💡 Note: a system PyMOL.app was detected, but GLINT is using the standalone conda PyMOL${NC}"
        echo -e "${BLUE}   This keeps dependency compatibility intact and does not affect your system PyMOL${NC}"
    fi
else
    echo -e "${YELLOW}⚠️  Using the system PyMOL.app (dependency conflicts are possible)${NC}"
    echo -e "${YELLOW}   If you encounter issues, reinstalling with conda PyMOL is recommended${NC}"
fi
echo ""
echo -e "${BLUE}📝 Note: the project has been renamed from GlueTK to GLINT${NC}"
echo ""
