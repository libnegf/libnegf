#!/bin/bash
# Author: Christoph Conrads (Forschungszentrum Jülich)
#
# The git pre-commit hook checks if Python and C++ code is properly formatted.
#
# The code requires YAPF for Python and clang-format for C/C++/CUDA.

set -e
set -u


# check for presence of formatting tools
if ! which clang-format >/dev/null
then
	>/dev/stderr echo 'clang-format not found. aborting.'
	exit 1
fi

if ! python3 -c 'import yapf' 2>/dev/null
then
	>/dev/stderr echo 'Python3 module YAPF not found. aborting.'
	exit 1
fi


needs_format=false


# check files
files="$(git diff --cached --name-only --diff-filter=ACM)"

while IFS='' read -r file
do
	if [[ "$file" =~ [.]([ch]pp|cu)$ ]]
	then
		if ! $(clang-format --dry-run --Werror -- "$file" 2>/dev/null)
		then
			echo "$file must be formatted with clang-format"
			needs_format=true
		fi
	elif [[ "$file" =~ [.]py$ ]]
	then
		if ! $(python3 -m yapf --quiet -- "$file")
		then
			echo "$file must be formatted with YAPF"
			needs_format=true
		fi
	fi
done <<<"$files"

if [ "$needs_format" = true ]
then
	exit 1
fi
