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

# Space separated folders in $ROOT_DIR
EXCLUDED_FOLDERS=".env .mypy_cache"

function find_cmd() {

  # Get path
  path=$1
  shift

  # Get modules
  modules=$(cat ${ROOT_DIR}/.gitmodules 2>/dev/null | grep path | sed 's/[[:space:]]*//g' | awk -F "=" -v root_dir=${ROOT_DIR} '{ print root_dir"/"$2 }' | xargs realpath | tr '\n' '|')
  if [[ ! -z $modules ]]; then
    modules=${modules::$((${#modules}-1))}
    modules="$modules|$(realpath ${ROOT_DIR}/.git)"
  fi

  # Get excluded folders
  excluded_folders=$(echo $EXCLUDED_FOLDERS | sed 's/ \+/|/g' | xargs realpath)

  if [[ ! -z $modules && ! -z $excluded_folders ]]; then
    find "$path" "$@" | xargs realpath | grep -Ev ^"$modules|$excluded_folders"
  elif [[ ! -z $modules ]]; then
    find "$path" "$@" | xargs realpath | grep -Ev ^$modules
  elif [[ ! -z $excluded_folders ]]; then
    find "$path" "$@" | xargs realpath | grep -Ev ^$excluded_folders
  else
    find "$path" "$@" | xargs realpath
  fi
}

# Make a list of files that need formatting
malformed_files=()
malformed_py_files=()

# For all files in this directory and all subdirectories...
for filename in $(find_cmd ${ROOT_DIR}/ -type f -iname "*.h" -or -iname "*.cpp"); do
  filename=$(realpath $(dirname $filename))/$(basename $filename)
  echo "Checking: $filename" >&2
  # ...check if there are any changes required.
  if ${CXX_CHECKER} --style=file --output-replacements-xml "$filename" | grep -q "<replacement "; then
    # This file requires changes, add it to the list.
    malformed_files=("$filename" ${malformed_files[@]})
  fi
done

for filename in $(find_cmd ${ROOT_DIR}/ -type f -iname "*.py"); do
  filename=$(realpath $(dirname $filename))/$(basename $filename)
  echo "Checking: $filename" >&2
  # ...check if there are any changes required.
  if [[ $(${PY_CHECKER} -d "$filename" | wc -l) > 0 ]]; then
    # This file requires changes, add it to the list.
    malformed_py_files=("$filename" ${malformed_py_files[@]})
  fi
done

# If any files require formatting, list them and return an error.
status=0

echo >&2
if ! [ ${#malformed_files[@]} -eq 0 ]; then
  echo "C++ files require formatting: ${malformed_files[@]}"    >&2
  echo                                                          >&2
  echo "Run the following command to auto-format these files:"  >&2
  echo "${CXX_CHECKER} --style=file -i ${malformed_files[@]}"   >&2
  echo                                                          >&2
  status=1
else
  echo "All C++ files are formatted correctly."                 >&2
  echo                                                          >&2
fi

if ! [ ${#malformed_py_files[@]} -eq 0 ]; then
  echo "Python files require formatting: ${malformed_py_files[@]}"  >&2
  echo                                                              >&2
  echo "Run the following command to auto-format these files:"      >&2
  echo "${PY_CHECKER} -i ${malformed_py_files[@]}"                  >&2
  echo                                                              >&2
  status=1
else
  echo "All Python files are formatted correctly."                  >&2
  echo                                                              >&2
fi

exit $status
