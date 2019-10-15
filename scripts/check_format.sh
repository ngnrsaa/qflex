#!/bin/bash

ROOT_DIR="$(realpath $(dirname $0))/../"

function find_cmd() {
  # Get path
  path=$1
  shift
  # Get modules
  modules=$(cat ${ROOT_DIR}/.gitmodules 2>/dev/null | grep path | sed 's/[[:space:]]*//g' | awk -F "=" -v root_dir=${ROOT_DIR} '{ print root_dir"/"$2 }' | tr '\n' '|')
  if [[ ! -z $modules ]]; then
    modules=${modules::-1}
    find "$path" "$@" | grep -Ev ^"$modules|${ROOT_DIR}/.git"
  else
    find "$path" "$@" 
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
  if clang-format --style=Google --output-replacements-xml "$filename" | grep -q "<replacement "; then
    # This file requires changes, add it to the list.
    malformed_files=("$filename" ${malformed_files[@]})
  fi
done

for filename in $(find_cmd ${ROOT_DIR}/ -type f -iname "*.py"); do
  filename=$(realpath $(dirname $filename))/$(basename $filename)
  echo "Checking: $filename" >&2
  # ...check if there are any changes required.
  if [[ $(yapf3 --style=Google -d "$filename" | wc -l) > 0 ]]; then
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
  echo "clang-format --style=Google -i ${malformed_files[@]}"   >&2
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
  echo "yapf3 --style=Google -i ${malformed_py_files[@]}"           >&2
  echo                                                              >&2
  status=1
else
  echo "All Python files are formatted correctly."                  >&2
  echo                                                              >&2
fi

exit $status
