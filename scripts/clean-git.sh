#!/bin/bash

# This script enables the "make clean" operation.

function get_path() {
  echo $(realpath $(dirname $1))
}

if [[ ! ($# == 0 || ($# == 1 && $1 == "-y")) ]]; then
  echo -e "\n\tUsage: $(get_path $0)/$(basename $0) [-y]\n" >&2
  echo -e "\tOptions:"                                      >&2
  echo -e "\t\t-y: Move ignored files if any\n\n"           >&2
  exit 1
fi

# Define root folder
ROOT_DIR=$(get_path $(get_path $0))

# Create temporary folder
TMP_DIR=/tmp/qflex-git_${RANDOM}${RANDOM}

# Get ignored files and move them to a temporary directory
IGNORED_FILES=$(git --work-tree=$ROOT_DIR status --ignored -s | grep ^'!!' | sed 's/\!\! //g' | xargs -I{} echo $ROOT_DIR/{})

if [[ -z $IGNORED_FILES ]]; then exit 0; fi

if [[ $1 != "-y" ]]; then
  echo "The following files are ignored by git:"                  >&2
  echo $IGNORED_FILES | tr ' ' '\n' | xargs -I{} echo -e "\t* {}" >&2
  echo -n "Move them to ($TMP_DIR)? (y/N) "                       >&2
  read ANSWER
else
  ANSWER="Y"
fi

if [[ $ANSWER == "Y" || $ANSWER == "y" ]]; then
  mkdir -vp $TMP_DIR >&2
  for file in $IGNORED_FILES; do
    file=$(get_path $file)/$(basename $file)
    dest=$(get_path $file)
    dest=${dest/$ROOT_DIR/$TMP_DIR}
    mkdir -vp $dest >&2
    mv -v $file $dest >&2
  done
  echo "Ignored files moved to: $TMP_DIR" >&2
  exit 0
else
  exit 1
fi
