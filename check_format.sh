#!/bin/bash
# Make a list of files that need formatting
malformed_files=()

# For all files in this directory and all subdirectories...
for filename in $(find . -iname "*.h" -or -iname "*.cpp"); do
  # ...check if there are any changes required.
  if clang-format --style=Google --output-replacements-xml "$filename" | grep -q "<replacement "; then
    # This file requires changes, add it to the list.
    malformed_files=("$filename" ${malformed_files[@]})
  fi
done

# If any files require formatting, list them and return an error.
if ! [ ${#malformed_files[@]} -eq 0 ]; then
  echo "Files require formatting: ${malformed_files[@]}"
  echo
  echo "Run the following command to auto-format these files:"
  echo "clang-format --style=Google -i ${malformed_files[@]}"
  echo
  exit 1
else
  echo "All files are formatted correctly."
  echo
  exit 0
fi
