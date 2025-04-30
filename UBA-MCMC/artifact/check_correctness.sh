#!/bin/bash

# Check if an argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <number-of-random-LMCs>"
    exit 1
fi

# Assign the first argument to a variable
num="$1"

# Check if the input is a valid positive integer (greater than 0)
if ! [[ "$num" =~ ^[0-9]+$ ]] || [ "$num" -le 0 ]; then
    echo "Error: Input must be a positive integer greater than 0."
    exit 1
fi

error=0 
valid_lmcs=0  # Counter for valid LMCs
lmcfile=lmc0.pm
file1="results/uba-test.txt"
file2="results/gfg-test.txt"

# Generate $num random LMCs
while [ "$valid_lmcs" -lt "$num" ]; do

  rm -f $file1 $file2 $lmcfile

  # Generate a random lmc
  python3 dtmc-gen.py $lmcfile

  echo "Check LMC $valid_lmcs..."

  # Prism-UBA model check
  ../prism/bin/prism -javamaxmem 200g -javastack 1g -timeout 30m -explicit -ltluba -ubaverbosity  1 -ubapure lmc0.pm -pf 'P=?[ HOA: {"'"UBAs/uba-3.hoa"'", "sigma" <- "sigma", "pi" <- "pi", "hash" <- "hash", "dollar" <- "dollar"}]' > $file1

  # Prism-GFG model check
  ../prism/bin/prism -javamaxmem 100g -javastack 1g -timeout 30m -explicit -gfgmc -ubaverbosity 1 lmc0.pm -pf 'P=?[ HOA: {"'"GFGs/gfg-3.hoa"'", "sigma_0" <- "sigma", "pi_0" <- "pi", "hash_0" <- "hash","dollar_0" <- "dollar" }]' > $file2

  # Read result line from file1
  line1=$(grep "^Result:" $file1)

  # Extract the value between "Result:" and "("
  value1=$(echo "$line1" | sed -E 's/^Result:[[:space:]]*([^()]*)\(.*$/\1/' | xargs)

  # Read result line from file2
  line2=$(grep "^Result:" $file2)

  # Extract the value between "Result:" and "("
  value2=$(echo "$line2" | sed -E 's/^Result:[[:space:]]*([^()]*)\(.*$/\1/' | xargs)

  # Check if sed produced any output
  if [ -z "$value2" ]; then
    echo "Parse error: Could not extract a value from line: '$line2'"
    continue  # Skip this iteration
  fi

  # Optionally: Check if it's a valid number (adjust as needed)
  if ! [[ "$value2" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Parse error: Extracted value is not a valid number: '$value2'"
    continue  # Skip this iteration
  fi

  echo "$value1, $value2"
  epsilon=0.00001

  diff=$(echo "$value1 - $value2" | bc -l)
  abs_diff=$(echo "$diff" | awk '{print ($1<0)?-$1:$1}')

  if (( $(echo "$abs_diff > $epsilon" | bc -l) )); then
    echo "Results for case n = $valid_lmcs do not match! Please check lmc0.pm."
    error=1
    break
  else
    # Increment the valid LMCs counter
    valid_lmcs=$((valid_lmcs + 1))
  fi

done

if [ "$error" != 0 ]; then
  echo "Error!"
else
  rm -f $file1 $file2 $lmcfile
  echo "All $num random LMCs match!"
fi
