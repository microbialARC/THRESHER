#!/usr/bin/env bash
set -e
# Parse command line arguments for forece reinstall option
FORCE_REINSTALL=false
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--force)
            FORCE_REINSTALL=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -f, --force    Force reinstall in existing thresher environment"
            echo "  -h, --help     Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

echo "Starting Thresher installation..."
# Detect if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Conda is not installed. Please install Miniconda or Anaconda first."
    exit 1
fi

# Initialize conda for bash shell
if command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)" || . "$(conda info --base)/etc/profile.d/conda.sh"
fi

# Get the directory of the install.sh script
INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Check if thresher environment exists
ENV_EXISTS=false
if conda env list | grep -q '^thresher\s'; then
    ENV_EXISTS=true
fi

# Handle installation based on environment existence and force flag
if [ "$ENV_EXISTS" = true ] && [ "$FORCE_REINSTALL" = false ]; then
    echo "Detected existing conda environment 'thresher' but --force flag not set."
    echo "Quitting installation to avoid overwriting existing environment."
    echo "If the environment is corrupted or dependencies are missing, please remove the existing thresher environment and rerun this script."
    echo "If you want to reinstall thresher by overwriting the existing environment, use:"
    echo "bash install.sh --force"
    exit 1
elif [ "$ENV_EXISTS" = true ] && [ "$FORCE_REINSTALL" = true ]; then
    echo "Force reinstalling in existing conda environment 'thresher'..."
elif [ "$ENV_EXISTS" = false ] && [ "$FORCE_REINSTALL" = true ]; then
    echo "--force flag detected but thresher environment does not exist."
    echo "Creating conda environment 'thresher'..."
    conda env create -f "$INSTALL_DIR/thresher.yml"
elif [ "$ENV_EXISTS" = false ] && [ "$FORCE_REINSTALL" = false ]; then
    # normal first install
    echo "Creating conda environment 'thresher'..."
    conda env create -f "$INSTALL_DIR/thresher.yml"
fi

# Activate the thresher conda environment
printf "\nActivating thresher environment...\n"
conda activate thresher

# Double check if the correct conda environment is activated
CURRENT_ENV=$(basename "$CONDA_PREFIX")
if [ "$CURRENT_ENV" != "thresher" ]; then
    echo "Wrong conda environment detected. Current environment: $CURRENT_ENV"
    echo "Please activate the thresher environment: conda activate thresher"
    exit 1
fi

# Install Thresher with verbose info
echo "Installing thresher package..."
pip install "$INSTALL_DIR" -v

# Print success message
cols=$(tput cols 2>/dev/null || echo "${COLUMNS:-80}")
printf '%*s\n' "$cols" '' | tr ' ' '='

if color="$(tput setaf 3 2>/dev/null)"; then
    reset="$(tput sgr0 2>/dev/null)"
else
    color=$'\033[0;33m'
    reset=$'\033[0m'
fi

printf '%b\n' "${color}Thresher installed successfully!${reset}"
printf '%b\n' "${color}Activate the environment with: conda activate thresher${reset}"
printf '%b\n' "${color}Then run 'thresher -h' to see the help page.${reset}"
printf '%*s\n' "$cols" '' | tr ' ' '='