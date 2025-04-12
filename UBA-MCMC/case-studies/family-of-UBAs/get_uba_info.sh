#!/bin/bash

# Check if a file argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <file>"
    exit 1
fi

file=$1

# Ensure the file exists
if [ ! -f "$file" ]; then
    echo "Error: File '$file' not found!"
    exit 1
fi

# Extract UBA number of states 
uba=$(grep "UBA has" "$file" | grep -oE '[0-9]+')

# Extract number of states in product
prod=$(awk '/product has/ { for (i=1; i<=NF; i++) if ($i == "has") print $(i+1) }' "$file")


# Extract the time in seconds from the line containing "Time for model checking:"
mctime=$(grep "Time for model checking:" "$file" | grep -oE '[0-9]+(\.[0-9]+)?')


# Output the total time
echo "#states = $uba, #product = $prod, time = $mctime s"


