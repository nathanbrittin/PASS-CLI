#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Local build helper — Windows release binary only.
#
# Full cross-platform distribution (Windows, Linux, macOS Intel, macOS ARM)
# is handled automatically by GitHub Actions when you push a version tag:
#
#   git tag v1.3.0 && git push origin v1.3.0
#
# See .github/workflows/release.yml for the full release pipeline.
# ---------------------------------------------------------------------------
set -euo pipefail

if [ ! -f "Cargo.toml" ]; then
    echo "Error: run this script from the project root."
    exit 1
fi

echo "Building release binary for the current platform..."
cargo build --release

echo "Copying to dist/..."
mkdir -p dist
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" || "$OSTYPE" == "win32" ]]; then
    cp target/release/pass-cli.exe dist/pass-cli.exe
    echo "✓ dist/pass-cli.exe updated"
else
    cp target/release/pass-cli dist/pass-cli
    echo "✓ dist/pass-cli updated"
fi
