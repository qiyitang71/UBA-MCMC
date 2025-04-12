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

# Extract GFG number of states 
gfg=$(grep "GFG has" "$file" | grep -oE '[0-9]+')

# Extract number of states in product
prod=$(awk '/product has/ { for (i=1; i<=NF; i++) if ($i == "has") print $(i+1) }' "$file")

# Extract the time in seconds from line 21
time_line21=$(sed -n '21p' "$file" | grep -oE '[0-9]+(\.[0-9]+)?')

# Extract the time in seconds from the line containing "Time for model checking:"
time_model_checking=$(grep "Time for model checking:" "$file" | grep -oE '[0-9]+(\.[0-9]+)?')

# Ensure both times are valid floating-point numbers
if [[ ! $time_line21 =~ ^[0-9]+(\.[0-9]+)?$ ]] || [[ ! $time_model_checking =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: Could not extract valid times from the file."
    exit 1
fi

# Calculate the total time using awk for floating-point addition
total_time=$(awk -v t1="$time_line21" -v t2="$time_model_checking" 'BEGIN { print t1 + t2 }')

# Output the total time
echo "#states = $gfg, #product = $prod, total_time = $total_time s: $time_line21 (trans), $time_model_checking (mc)"


