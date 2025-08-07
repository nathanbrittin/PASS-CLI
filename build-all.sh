#!/usr/bin/env bash
set -euo pipefail

# Set environment variables for pkg-config
export PKG_CONFIG_ALLOW_SYSTEM_LIBS=1
export PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Project name (change this to match your actual binary name)
PROJECT_NAME="pass-cli"

echo -e "${BLUE}=== Cross-compilation build script for WSL ===${NC}"

# Check if we're in a Rust project
if [ ! -f "Cargo.toml" ]; then
    echo -e "${RED}Error: Cargo.toml not found. Are you in a Rust project directory?${NC}"
    exit 1
fi

# 1) Install required dependencies for cross-compilation in WSL
echo -e "${YELLOW}Checking WSL dependencies...${NC}"

# Check if we need to install dependencies
NEED_INSTALL=false

if ! command -v x86_64-w64-mingw32-gcc &> /dev/null; then
    NEED_INSTALL=true
fi

if ! command -v cmake &> /dev/null; then
    NEED_INSTALL=true
fi

if [ "$NEED_INSTALL" = true ]; then
    echo -e "${YELLOW}Installing cross-compilation dependencies...${NC}"
    sudo apt update
    sudo apt install -y \
        gcc-mingw-w64-x86-64 \
        gcc-mingw-w64-i686 \
        g++-mingw-w64-x86-64 \
        g++-mingw-w64-i686 \
        cmake \
        build-essential \
        pkg-config \
        libssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev
    
    echo -e "${GREEN}✓ Dependencies installed${NC}"
else
    echo -e "${GREEN}✓ All dependencies already installed${NC}"
fi

# 2) Ensure Rust knows about our targets
echo -e "${YELLOW}Adding Rust targets...${NC}"
rustup target add x86_64-unknown-linux-gnu
rustup target add x86_64-pc-windows-gnu

# 3) Set up cross-compilation environment for Windows
export CC_x86_64_pc_windows_gnu=x86_64-w64-mingw32-gcc
export CXX_x86_64_pc_windows_gnu=x86_64-w64-mingw32-g++
export AR_x86_64_pc_windows_gnu=x86_64-w64-mingw32-ar

# 4) Build for Windows (GNU): ${PROJECT_NAME}.exe
echo -e "${BLUE}Building Windows binary...${NC}"
if cargo build --release --target x86_64-pc-windows-gnu; then
    echo -e "${GREEN}✓ Windows build successful${NC}"
else
    echo -e "${RED}✗ Windows build failed${NC}"
    exit 1
fi

# 5) Build for Linux (glibc): ${PROJECT_NAME}
echo -e "${BLUE}Building Linux binary...${NC}"
if cargo build --release --target x86_64-unknown-linux-gnu; then
    echo -e "${GREEN}✓ Linux build successful${NC}"
else
    echo -e "${RED}✗ Linux build failed${NC}"
    exit 1
fi

# 6) Package outputs
echo -e "${BLUE}Packaging binaries...${NC}"
rm -rf dist && mkdir -p dist

# Copy Windows binary
if [ -f "target/x86_64-pc-windows-gnu/release/${PROJECT_NAME}.exe" ]; then
    cp "target/x86_64-pc-windows-gnu/release/${PROJECT_NAME}.exe" dist/
    echo -e "${GREEN}✓ Windows binary packaged${NC}"
else
    echo -e "${RED}✗ Windows binary not found${NC}"
    exit 1
fi

# Copy Linux binary
if [ -f "target/x86_64-unknown-linux-gnu/release/${PROJECT_NAME}" ]; then
    cp "target/x86_64-unknown-linux-gnu/release/${PROJECT_NAME}" dist/
    echo -e "${GREEN}✓ Linux binary packaged${NC}"
else
    echo -e "${RED}✗ Linux binary not found${NC}"
    exit 1
fi

# 7) Make Linux binary executable (just in case)
chmod +x "dist/${PROJECT_NAME}"

# 8) Create release archives
echo -e "${BLUE}Creating release archives...${NC}"

# Get version from Cargo.toml if available
VERSION=""
if command -v grep &> /dev/null && [ -f "Cargo.toml" ]; then
    VERSION=$(grep '^version = ' Cargo.toml | head -n1 | cut -d'"' -f2)
    if [ -n "$VERSION" ]; then
        VERSION="v${VERSION}"
    fi
fi

# Fallback to git tag or timestamp
if [ -z "$VERSION" ]; then
    if command -v git &> /dev/null && git rev-parse --git-dir > /dev/null 2>&1; then
        VERSION=$(git describe --tags --abbrev=0 2>/dev/null || echo "")
    fi
fi

if [ -z "$VERSION" ]; then
    VERSION=$(date +"%Y%m%d-%H%M%S")
fi

# Create archives directory
mkdir -p releases

# Create Windows ZIP archive
WIN_ARCHIVE="${PROJECT_NAME}-${VERSION}-windows-x64.zip"
echo -e "${YELLOW}Creating Windows archive: ${WIN_ARCHIVE}${NC}"

# Create temporary directory for Windows release
WIN_TEMP_DIR="temp_win_release"
mkdir -p "$WIN_TEMP_DIR"
cp "dist/${PROJECT_NAME}.exe" "$WIN_TEMP_DIR/"

# Add README if it exists
if [ -f "README.md" ]; then
    cp "README.md" "$WIN_TEMP_DIR/"
fi

# Add LICENSE if it exists
if [ -f "LICENSE" ] || [ -f "LICENSE.txt" ] || [ -f "LICENSE.md" ]; then
    find . -maxdepth 1 -name "LICENSE*" -exec cp {} "$WIN_TEMP_DIR/" \;
fi

# Create ZIP (using zip command if available, otherwise use tar)
if command -v zip &> /dev/null; then
    (cd "$WIN_TEMP_DIR" && zip -r "../releases/$WIN_ARCHIVE" .)
else
    echo -e "${YELLOW}zip command not found, installing...${NC}"
    sudo apt install -y zip
    (cd "$WIN_TEMP_DIR" && zip -r "../releases/$WIN_ARCHIVE" .)
fi

# Clean up temp directory
rm -rf "$WIN_TEMP_DIR"

if [ -f "releases/$WIN_ARCHIVE" ]; then
    echo -e "${GREEN}✓ Windows archive created: releases/${WIN_ARCHIVE}${NC}"
else
    echo -e "${RED}✗ Failed to create Windows archive${NC}"
fi

# Create Linux TAR.GZ archive
LINUX_ARCHIVE="${PROJECT_NAME}-${VERSION}-linux-x64.tar.gz"
echo -e "${YELLOW}Creating Linux archive: ${LINUX_ARCHIVE}${NC}"

# Create temporary directory for Linux release
LINUX_TEMP_DIR="temp_linux_release"
mkdir -p "$LINUX_TEMP_DIR"
cp "dist/${PROJECT_NAME}" "$LINUX_TEMP_DIR/"

# Add README if it exists
if [ -f "README.md" ]; then
    cp "README.md" "$LINUX_TEMP_DIR/"
fi

# Add LICENSE if it exists
if [ -f "LICENSE" ] || [ -f "LICENSE.txt" ] || [ -f "LICENSE.md" ]; then
    find . -maxdepth 1 -name "LICENSE*" -exec cp {} "$LINUX_TEMP_DIR/" \;
fi

# Create TAR.GZ
tar -czf "releases/$LINUX_ARCHIVE" -C "$LINUX_TEMP_DIR" .

# Clean up temp directory
rm -rf "$LINUX_TEMP_DIR"

if [ -f "releases/$LINUX_ARCHIVE" ]; then
    echo -e "${GREEN}✓ Linux archive created: releases/${LINUX_ARCHIVE}${NC}"
else
    echo -e "${RED}✗ Failed to create Linux archive${NC}"
fi

# 9) Display results
echo -e "${GREEN}=== Build Complete! ===${NC}"
echo -e "${BLUE}Binaries in ./dist:${NC}"
ls -lh dist/

echo -e "\n${BLUE}Release archives in ./releases:${NC}"
ls -lh releases/

# Show archive contents for verification
echo -e "\n${BLUE}Archive contents:${NC}"
if [ -f "releases/$WIN_ARCHIVE" ]; then
    echo -e "${YELLOW}Windows archive contents:${NC}"
    if command -v unzip &> /dev/null; then
        unzip -l "releases/$WIN_ARCHIVE" | tail -n +4 | head -n -2
    else
        echo "  ${PROJECT_NAME}.exe"
        [ -f "README.md" ] && echo "  README.md"
        find . -maxdepth 1 -name "LICENSE*" -printf "  %f\n"
    fi
fi

if [ -f "releases/$LINUX_ARCHIVE" ]; then
    echo -e "${YELLOW}Linux archive contents:${NC}"
    tar -tzf "releases/$LINUX_ARCHIVE" | sed 's/^/  /'
fi

echo -e "\n${GREEN}Cross-compilation and packaging successful! 🎉${NC}"
echo -e "${BLUE}Ready for distribution!${NC}"