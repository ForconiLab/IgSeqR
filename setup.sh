#!/bin/bash
VERSION="v0.1.0"

# Function to print messages
print_message() {
    echo -e "\n$1\n"
}

# Function to print commands
print_command() {
    echo -ne "\033[1;34m$1...\033[0m"  # Blue color for commands
}

# Function to print errors
print_error() {
    echo -e "\033[1;31m$1\033[0m"  # Red color for errors
}

# Function to print interaction prompts
print_prompt() {
    echo -e "\033[1;33m$1\033[0m"  # Yellow color for interaction prompts
}

echo ==============================================================
echo "IgSeqR $VERSION Installer"
echo ==============================================================

# Confirm installation
print_prompt "This will install IgSeqR and its dependencies. Do you want to continue? (y/n): "
read choice
case "$choice" in 
  y|Y ) print_command "Proceeding with installation..."; echo " DONE";;
  n|N ) print_error "Installation aborted."; exit 1;;
  * ) print_error "Invalid input. Installation aborted."; exit 1;;
esac

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    print_error "Conda is not installed. Please install Conda and rerun the setup script."
    exit 1
fi
print_command "Installing Conda environment"
# Check if the Conda environment already exists
ENV_NAME="igseqr"
if conda env list | grep -q "$ENV_NAME"; then
    print_error " A Conda environment named '$ENV_NAME' already exists."
    print_prompt "Do you want to remove the existing environment and create a new one? (y/n): "
    read remove_choice
    case "$remove_choice" in 
      y|Y ) 
        print_command "Removing existing Conda environment"
        conda env remove -n "$ENV_NAME"
        if [ $? -eq 0 ]; then
            echo " DONE"
        else
            print_error "Failed to remove existing Conda environment."
            exit 1
        fi
        ;;
      n|N ) 
        print_prompt "Enter a new name for the Conda environment: "
        read new_env_name
        ENV_NAME="$new_env_name"
        ;;
      * ) 
        print_error "Invalid input. Installation aborted."
        exit 1
        ;;
    esac
fi

# Install Conda environment
conda env create --name "$ENV_NAME" --file=environment.yml
if [ $? -eq 0 ]; then
    echo " DONE"
else
    print_error "Failed to install Conda environment. Please check the environment.yml file and try again."
    exit 1
fi

# Activate the Conda environment
print_command "Activating Conda environment"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$ENV_NAME"
if [ $? -eq 0 ]; then
    echo " DONE"
else
    print_error "Failed to activate Conda environment. Please check your Conda installation."
    exit 1
fi

# Determine the Conda environment directory
CONDA_ENV_DIR=$(conda info --base)/envs/"$ENV_NAME"
print_message "Using Conda environment directory: $CONDA_ENV_DIR"

# Create necessary directories inside the Conda environment
print_command "Creating necessary directories inside Conda environment"
mkdir -p "$CONDA_ENV_DIR/bin/data/igseqr"
if [ $? -eq 0 ]; then
    echo " DONE"
else
    print_error "Failed to create necessary directories."
    exit 1
fi

# Copy the igseqr.sh script to the Conda environment bin directory
if [ -f "$(pwd)/bin/igseqr.sh" ]; then
    print_command "Copying igseqr.sh script to Conda environment bin directory"
    cp $(pwd)/bin/igseqr.sh "$CONDA_ENV_DIR/bin/"
    chmod +x "$CONDA_ENV_DIR/bin/igseqr.sh"
    if [ $? -eq 0 ]; then
        echo " DONE"
    else
        print_error "Failed to copy igseqr.sh script."
        exit 1
    fi
else
    print_error "igseqr.sh script not found. Please check your repository structure."
    exit 1
fi

# Link the igseqr.sh script to the Conda environment bin directory
print_command "Linking igseqr.sh script to Conda environment bin directory"
ln -sf "$CONDA_ENV_DIR/bin/igseqr.sh" "$CONDA_ENV_DIR/bin/igseqr"
if [ $? -eq 0 ]; then
    echo " DONE"
else
    print_error "Failed to link igseqr.sh script."
    exit 1
fi

# Copy data files to the Conda environment bin/data directory
if [ -d "$(pwd)/data" ]; then
    print_command "Copying data files to Conda environment bin/data directory"
    cp -r $(pwd)/data/* "$CONDA_ENV_DIR/bin/data/igseqr"
    if [ $? -eq 0 ]; then
        echo " DONE"
    else
        print_error "Failed to copy data files."
        exit 1
    fi
else
    print_error "Data directory not found. Please check your repository structure."
    exit 1
fi

# Check if all necessary directories and files are in place
print_command "Verifying file structure..."
if [ -f "$CONDA_ENV_DIR/bin/igseqr" ] && [ -d "$CONDA_ENV_DIR/bin/data/igseqr" ]; then
    echo " DONE"
else
    print_error "Some directories or files are missing. Please check your repository structure."
    exit 1
fi

print_message "Installation complete! You can now run 'igseqr' from anywhere when the '$ENV_NAME' Conda environment is active."
