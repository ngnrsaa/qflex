#!/bin/bash

# This script checks the formatting of all C++ and python files in the
# repository. If issues are found, a command to automatically resolve them will
# be printed as output.

RESET="\e[0m"
RED="\e[91m"
GREEN="\e[92m"
CYAN="\e[96m"

ROOT_DIR="$(realpath $(realpath $(dirname $0))/../)"

PY_CHECKER=yapf
PY_CHECKER_VERSION=$(cat $ROOT_DIR/scripts/requirements.txt | grep ^${PY_CHECKER}== | cut -d '=' -f 3)

CXX_CHECKER=clang-format
CXX_CHECKER_VERSION=$(cat $ROOT_DIR/scripts/requirements.txt | grep ^${CXX_CHECKER}== | cut -d '=' -f 3)

if ! which $PY_CHECKER >/dev/null || [[ $($PY_CHECKER --version | cut -d ' ' -f 2) != $PY_CHECKER_VERSION ]]; then
  echo -e "$PY_CHECKER is not available. Install $PY_CHECKER as:\n$ python3 -m pip install $PY_CHECKER==${PY_CHECKER_VERSION}" >&2
  exit 1
fi

if ! which $CXX_CHECKER >/dev/null || [[ $($CXX_CHECKER --version | cut -d ' ' -f 3) != $CXX_CHECKER_VERSION ]]; then
  echo -e "$CXX_CHECKER is not available. Install $CXX_CHECKER as:\n$ python3 -m pip install $CXX_CHECKER==${CXX_CHECKER_VERSION}" >&2
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
    echo -ne "${CYAN}[    ] Checking: $filename${RESET}" >&2
    # ...check if there are any changes required.
    if ${CXX_CHECKER} --style=google --output-replacements-xml "$filename" | grep -q "<replacement "; then
      # This file requires changes, add it to the list.
      echo -ne '"'$filename'" '
      echo -ne "\r${CYAN}[${RED}FAIL${RESET}" >&2
    else
      echo -ne "\r${CYAN}[${GREEN}PASS${RESET}" >&2
    fi
    echo >&2
  done
}

function check_py_format {
  while read filename; do
    filename=$(realpath $(dirname "$filename"))/$(basename "$filename")
    echo -ne "${CYAN}[    ] Checking: $filename${RESET}" >&2
    # ...check if there are any changes required.
    if [[ $(${PY_CHECKER} --style=google -d "$filename" | wc -l) > 0 ]]; then
      # This file requires changes, add it to the list.
      echo -ne '"'$filename'" '
      echo -ne "\r${CYAN}[${RED}FAIL${RESET}" >&2
    else
      echo -ne "\r${CYAN}[${GREEN}PASS${RESET}" >&2
    fi
    echo >&2
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
  echo "${CXX_CHECKER} --style=google -i ${malformed_files}"   >&2
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
  echo "${PY_CHECKER} --style=google -i ${malformed_py_files}"   >&2
  echo                                                           >&2
  status=1
else
  echo "All Python files are formatted correctly."               >&2
  echo                                                           >&2
fi

exit $status
