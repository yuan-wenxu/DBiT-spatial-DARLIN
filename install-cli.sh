#!/bin/bash
set -euo pipefail

show_help() {
    cat <<'EOF'
Usage: ./install-cli.sh

Install the dbit command for the current user by creating a symbolic link in
~/.local/bin. Set DBIT_INSTALL_DIR to use a different destination directory.

If needed, this script adds ~/.local/bin to PATH in ~/.bashrc.
EOF
}

if [[ $# -gt 1 ]]; then
    echo "Error: install-cli.sh accepts no positional arguments." >&2
    exit 1
fi
if [[ $# -eq 1 ]]; then
    case "$1" in
        -h|--help) show_help; exit 0 ;;
        *) echo "Error: unknown option or argument '$1'." >&2; exit 1 ;;
    esac
fi

if [[ -z ${HOME:-} && -z ${DBIT_INSTALL_DIR:-} ]]; then
    echo "Error: HOME is not set; set DBIT_INSTALL_DIR explicitly." >&2
    exit 1
fi

INSTALLER_PATH=$(readlink -f "${BASH_SOURCE[0]}")
REPO_DIR=$(cd "$(dirname "$INSTALLER_PATH")" && pwd)
SOURCE="$REPO_DIR/script/dbit.sh"
INSTALL_DIR=${DBIT_INSTALL_DIR:-$HOME/.local/bin}

mkdir -p "$INSTALL_DIR"
INSTALL_DIR=$(cd "$INSTALL_DIR" && pwd)
TARGET="$INSTALL_DIR/dbit"

if [[ -e "$TARGET" || -L "$TARGET" ]]; then
    existing_target=$(readlink -f "$TARGET" 2>/dev/null || true)
    source_target=$(readlink -f "$SOURCE")
    if [[ "$existing_target" == "$source_target" ]]; then
        echo "dbit is already installed at $TARGET"
    else
        echo "Error: $TARGET already exists and points elsewhere." >&2
        echo "Remove or rename it before installing dbit." >&2
        exit 1
    fi
else
    ln -s "$SOURCE" "$TARGET"
    echo "Installed dbit at $TARGET"
fi

bashrc_has_local_bin=false
if [[ -n ${HOME:-} && "$INSTALL_DIR" == "$(readlink -m "$HOME/.local/bin")" ]]; then
    BASHRC=${DBIT_BASHRC:-$HOME/.bashrc}
    if [[ -f "$BASHRC" ]] && awk -v home="$HOME" '
        /^[[:space:]]*#/ { next }
        index($0, "$HOME/.local/bin") ||
        index($0, "${HOME}/.local/bin") ||
        index($0, "~/.local/bin") ||
        index($0, home "/.local/bin") { found = 1 }
        END { exit !found }
    ' "$BASHRC"; then
        bashrc_has_local_bin=true
        echo "Found ~/.local/bin in $BASHRC"
    else
        echo "No ~/.local/bin PATH entry found in $BASHRC"
        {
            echo
            echo "# Added by DBiT-spatial-DARLIN install-cli"
            echo 'export PATH="$HOME/.local/bin:$PATH"'
        } >> "$BASHRC"
        bashrc_has_local_bin=true
        echo "Added ~/.local/bin to PATH in $BASHRC"
    fi
fi

if [[ ":$PATH:" == *":$INSTALL_DIR:"* ]]; then
    echo "Run 'dbit --help' to verify the installation."
elif $bashrc_has_local_bin; then
    echo "Reload the configured PATH with: source \"$BASHRC\""
else
    echo
    echo "Add this line to your shell startup file, then start a new shell:"
    printf 'export PATH=%q:$PATH\n' "$INSTALL_DIR"
fi
