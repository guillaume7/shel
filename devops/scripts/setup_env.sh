#!/bin/bash
# Script to set up a Python virtual environment for SHEL

# Configuration
ENV_NAME="venv"
PYTHON_VERSION="python3"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Print header
echo -e "${GREEN}=== Setting up SHEL Python environment ===${NC}"

# Check if Python 3 is installed
if ! command -v $PYTHON_VERSION &> /dev/null; then
    echo -e "${RED}Error: $PYTHON_VERSION is not installed. Please install Python 3 first.${NC}"
    exit 1
fi

# Check for virtual environment
if [ -d "$ENV_NAME" ]; then
    echo -e "${YELLOW}Warning: $ENV_NAME directory already exists.${NC}"
    read -p "Do you want to remove it and create a new one? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        rm -rf "$ENV_NAME"
    else
        echo "Using existing environment."
        source "$ENV_NAME/bin/activate"
        echo -e "${GREEN}Environment activated. Run 'deactivate' to exit.${NC}"
        exit 0
    fi
fi

# Create virtual environment
echo "Creating virtual environment..."
$PYTHON_VERSION -m venv "$ENV_NAME"

# Activate virtual environment
echo "Activating virtual environment..."
source "$ENV_NAME/bin/activate"

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install dependencies
echo "Installing dependencies..."
pip install -r requirements.txt

# Install PyQt5 separately to ensure it works
echo "Installing PyQt5..."
pip install PyQt5

# Run code generation script
echo "Running __init__.py file generation script..."
bash create_init_files.sh

# Install the package in development mode
echo "Installing SHEL in development mode..."
pip install -e .

echo -e "${GREEN}=== Setup complete! ===${NC}"
echo -e "Virtual environment created and activated."
echo -e "Run 'source $ENV_NAME/bin/activate' to activate the environment in the future."
echo -e "Run 'deactivate' to exit the virtual environment."
