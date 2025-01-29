#!/usr/bin/env bash
BIN_DIR="$CONDA_PREFIX/bin"
cp -r thresher/workflow/ $CONDA_PREFIX
cp thresher/bin/cli.py "$BIN_DIR/thresher"

chmod +x "$BIN_DIR/thresher"

if [[ ":$PATH:" != *":$BIN_DIR:"* ]]; then
    echo "export PATH=\"$BIN_DIR:\$PATH\"" >> "$HOME/.bashrc"
    echo "Added $BIN_DIR to PATH. Please restart your shell or run 'source ~/.bashrc'"
fi

echo "Thresher installed successfully! Run 'thresher -h' to start with help page."