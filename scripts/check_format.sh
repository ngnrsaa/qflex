#!/bin/bash

# This script checks the formatting of all C++ and python files in the
# repository. If issues are found, a command to automatically resolve them will
# be printed as output.

ROOT_DIR="$(realpath $(dirname $0))/../"
if which yapf3 > /dev/null; then
  PY_CHECKER=yapf3
elif which yapf > /dev/null; then
  PY_CHECKER=yapf
else
  echo "No available python checkers installed." >&2
  exit 1
fi
if which clang-format > /dev/null; then
  CXX_CHECKER=clang-format
else
  echo "No available c++ checkers installed." >&2
  exit 1
fi

# Find files within git repository given extension
function git_find() {
  while [[ $# > 0 ]]; do

    # Get extension
    ext=$1; shift

    # Get files
    git --git-dir=$ROOT_DIR/.git --work-tree=$ROOT_DIR ls-files --exclude-per-directory=.gitignore -co ${ROOT_DIR} | \
      awk -v root=$ROOT_DIR -v ext=$ext -F. '$NF == ext { print root"/"$0 }'

  done
}

function check_cxx_format {
  while read filename; do
    filename=$(realpath $(dirname "$filename"))/$(basename "$filename")
    echo "Checking: $filename" >&2
    # ...check if there are any changes required.
    if ${CXX_CHECKER} --style=file --output-replacements-xml "$filename" | grep -q "<replacement "; then
      # This file requires changes, add it to the list.
      echo -ne '"'$filename'" '
    fi
  done
}

function check_py_format {
  while read filename; do
    filename=$(realpath $(dirname "$filename"))/$(basename "$filename")
    echo "Checking: $filename" >&2
    # ...check if there are any changes required.
    if [[ $(${PY_CHECKER} -d "$filename" | wc -l) > 0 ]]; then
      # This file requires changes, add it to the list.
      echo -ne '"'$filename'" '
    fi
  done
}

malformed_files=$(git_find cpp h | check_cxx_format)
malformed_py_files=$(git_find py | check_py_format)

# If any files require formatting, list them and return an error.
status=0

echo >&2
if [[ -n ${malformed_files} ]]; then
  echo "C++ files require formatting: ${malformed_files}"      >&2
  echo                                                         >&2
  echo "Run the following command to auto-format these files:" >&2
  echo "${CXX_CHECKER} --style=file -i ${malformed_files}"     >&2
  echo                                                         >&2
  status=1
else
  echo "All C++ files are formatted correctly."                >&2
  echo                                                         >&2
fi

if [[ -n ${malformed_py_files} ]]; then
  echo "Python files require formatting: ${malformed_py_files}"  >&2
  echo                                                           >&2
  echo "Run the following command to auto-format these files:"   >&2
  echo "${PY_CHECKER} -i ${malformed_py_files}"                  >&2
  echo                                                           >&2
  status=1
else
  echo "All Python files are formatted correctly."               >&2
  echo                                                           >&2
fi

exit $status
