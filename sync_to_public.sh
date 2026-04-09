#!/bin/bash
# ============================================================================
# GLINT 私有仓库 → 公开仓库 同步脚本
# 用途：将开发版代码同步到开源仓库，排除内部文件
# ============================================================================

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 默认配置
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PRIVATE_REPO="${PRIVATE_REPO:-$SCRIPT_DIR}"
PUBLIC_REPO="${PUBLIC_REPO:-$(dirname "$PRIVATE_REPO")/GLINT}"
SYNC_EXCLUDE="${PRIVATE_REPO}/.sync_exclude"

# 帮助信息
show_help() {
    echo "用法: $0 [选项]"
    echo ""
    echo "选项:"
    echo "  -p, --public-repo DIR   公开仓库路径 (默认: ../GLINT)"
    echo "  -r, --private-repo DIR   私有仓库路径 (默认: 脚本所在目录)"
    echo "  -d, --dry-run            仅预览同步内容，不执行实际同步"
    echo "  -n, --no-commit          同步但不自动提交"
    echo "  -h, --help               显示此帮助信息"
    echo ""
    echo "环境变量:"
    echo "  PRIVATE_REPO             私有仓库路径"
    echo "  PUBLIC_REPO              公开仓库路径"
    echo ""
    echo "示例:"
    echo "  $0                                    # 使用默认路径同步"
    echo "  $0 -p ~/repos/GLINT                   # 指定公开仓库路径"
    echo "  $0 --dry-run                          # 预览同步内容"
}

# 参数解析
DRY_RUN=false
AUTO_COMMIT=true

while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--public-repo)
            PUBLIC_REPO="$2"
            shift 2
            ;;
        -r|--private-repo)
            PRIVATE_REPO="$2"
            shift 2
            ;;
        -d|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -n|--no-commit)
            AUTO_COMMIT=false
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo -e "${RED}错误: 未知选项 $1${NC}"
            show_help
            exit 1
            ;;
    esac
done

# 验证路径
if [[ ! -d "$PRIVATE_REPO" ]]; then
    echo -e "${RED}错误: 私有仓库路径不存在: $PRIVATE_REPO${NC}"
    exit 1
fi

if [[ ! -f "$SYNC_EXCLUDE" ]]; then
    echo -e "${RED}错误: 找不到排除配置文件: $SYNC_EXCLUDE${NC}"
    exit 1
fi

# 检查公开仓库是否存在
if [[ ! -d "$PUBLIC_REPO" ]]; then
    echo -e "${YELLOW}公开仓库不存在，正在创建: $PUBLIC_REPO${NC}"
    mkdir -p "$PUBLIC_REPO"
    cd "$PUBLIC_REPO"
    git init
    echo -e "${GREEN}✓ 已初始化公开仓库${NC}"
fi

# 检查公开仓库是否为 Git 仓库
if [[ ! -d "$PUBLIC_REPO/.git" ]]; then
    echo -e "${RED}错误: 公开仓库不是有效的 Git 仓库: $PUBLIC_REPO${NC}"
    exit 1
fi

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}GLINT 同步工具${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "私有仓库: ${GREEN}$PRIVATE_REPO${NC}"
echo -e "公开仓库: ${GREEN}$PUBLIC_REPO${NC}"
echo -e "排除配置: ${GREEN}$SYNC_EXCLUDE${NC}"
echo ""

# 预览模式
if [[ "$DRY_RUN" == true ]]; then
    echo -e "${YELLOW}[DRY-RUN] 预览将要同步的文件:${NC}"
    echo ""
    rsync -av --dry-run --exclude-from="$SYNC_EXCLUDE" \
        --delete --delete-excluded \
        "$PRIVATE_REPO/" "$PUBLIC_REPO/" | grep -v '/$' | head -50
    echo ""
    echo -e "${YELLOW}提示: 使用不带 --dry-run 参数执行实际同步${NC}"
    exit 0
fi

# 执行同步
echo -e "${BLUE}开始同步...${NC}"
rsync -av --exclude-from="$SYNC_EXCLUDE" \
    --delete --delete-excluded \
    "$PRIVATE_REPO/" "$PUBLIC_REPO/"

echo ""
echo -e "${GREEN}✓ 同步完成${NC}"

# 自动提交
if [[ "$AUTO_COMMIT" == true ]]; then
    echo ""
    echo -e "${BLUE}正在提交到公开仓库...${NC}"
    
    cd "$PUBLIC_REPO"
    
    # 检查是否有变更
    if [[ -z $(git status --porcelain) ]]; then
        echo -e "${YELLOW}没有变更需要提交${NC}"
        exit 0
    fi
    
    # 获取版本信息
    VERSION_FILE="$PRIVATE_REPO/glint/_version.py"
    if [[ -f "$VERSION_FILE" ]]; then
        VERSION=$(grep -oP '__version__\s*=\s*["\x27]\K[^"\x27]+' "$VERSION_FILE" 2>/dev/null || echo "unknown")
        COMMIT_MSG="Sync from GLINT_dev v$VERSION ($(date +%Y-%m-%d))"
    else
        COMMIT_MSG="Sync from GLINT_dev ($(date +%Y-%m-%d))"
    fi
    
    git add -A
    git commit -m "$COMMIT_MSG"
    
    echo -e "${GREEN}✓ 已提交: $COMMIT_MSG${NC}"
    
    # 提示推送
    echo ""
    echo -e "${YELLOW}提示: 运行以下命令推送到远程仓库:${NC}"
    echo -e "  cd $PUBLIC_REPO && git push origin main"
fi

echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${GREEN}同步完成！${NC}"
echo -e "${BLUE}========================================${NC}"
