#!/bin/bash

for file in $(find $1 -name "results.csv"); do
    if [ -f "$(dirname $file)/results.png" ]; then
        continue
    fi
    python3 "$(dirname $0)/results.py" $file -o "$(dirname $file )/results.png"
done