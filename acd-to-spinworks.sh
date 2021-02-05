#!/usr/bin/env bash

for file in "$@"; do
    grep -E '^Frequency|^[0-9-]' < "$file" | sed 's/Frequency = /ti: '"\n"'##frequency /g' | sed 's/\t/  /' | sed 's/,/./g' | tr -d '\015' > "$file.spin"
done
