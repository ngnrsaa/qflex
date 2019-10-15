#!/bin/bash

# Define root folder
ROOT_DIR="$(realpath $(dirname $0))/../"

# Create temporary folder
TMP_DIR=/tmp/qflex-git_${RANDOM}${RANDOM}
mkdir -p $TMP_DIR

# Get ignored files and move them to a temporary directory
git --work-tree=$ROOT_DIR status --ignored -s | grep ^'!!' | sed 's/\!\! //g' | xargs -I{} mv -v $ROOT_DIR/{} $TMP_DIR/

# Print temporary folder
echo "Ignored files are moved to: $TMP_DIR" >&2
