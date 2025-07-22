#!/usr/bin/env bash
set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Project name (change this to match your actual binary name)
PROJECT_NAME="PASS-CLI"

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
        libssl-dev
    
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

# 8) Display results
echo -e "${GREEN}=== Build Complete! ===${NC}"
echo -e "${BLUE}Binaries ready in ./dist:${NC}"
ls -lh dist/

# Optional: Show file info
echo -e "\n${BLUE}File details:${NC}"
file dist/*

# Optional: Show sizes
echo -e "\n${BLUE}Binary sizes:${NC}"
du -h dist/*

echo -e "\n${GREEN}Cross-compilation successful! 🎉${NC}"./buid