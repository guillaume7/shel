#!/bin/bash
# Create __init__.py files in all subdirectories with proper docstrings

BASE_DIR="/home/guillaume-riflet/git/shel/src/python"

# Create empty __init__.py files in all subdirectories
for dir in $(find $BASE_DIR -type d | grep -v "__pycache__" | sort); do
    # Skip the base directory
    if [ "$dir" != "$BASE_DIR" ]; then
        init_file="$dir/__init__.py"
        
        # Get directory name for docstring
        dir_name=$(basename "$dir")
        
        # Create __init__.py with appropriate docstring if it doesn't exist
        if [ ! -f "$init_file" ]; then
            echo "Creating $init_file"
            echo '"""' > "$init_file"
            echo "$dir_name module for SHEL." >> "$init_file"
            echo '"""' >> "$init_file"
        else
            echo "File already exists: $init_file"
        fi
    fi
done

# Create __init__.py in tests directory
tests_dir="/home/guillaume-riflet/git/shel/tests/python"
tests_init="$tests_dir/__init__.py"
if [ ! -f "$tests_init" ]; then
    echo "Creating $tests_init"
    echo '"""' > "$tests_init"
    echo "Tests for the SHEL Python implementation." >> "$tests_init"
    echo '"""' >> "$tests_init"
else
    echo "File already exists: $tests_init"
fi

echo "Done! All directories now have __init__.py files."
