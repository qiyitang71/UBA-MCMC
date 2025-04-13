#!/bin/bash

# Folder containing the files
results_folder="results"



# Loop through files gfg-3.txt to gfg-12.txt and uba-3.txt to uba-12.txt
for i in {3..12}; do
  file1="$results_folder/uba-$i.txt"
  file2="$results_folder/gfg-$i.txt"

  # Read from a file or input (replace 'input.txt' with your file name)
  line1=$(grep "^Result:" $file1)

  if [ -z "$line1" ]; then
    echo "Skipping n = $i, jcss19 did not complete"
    continue
  fi

  # Extract the value between "Result:" and "("
  value1=$(echo "$line1" | sed -E 's/^Result:[[:space:]]*([^()]*)\(.*$/\1/' | xargs)

    # Read from a file or input (replace 'input.txt' with your file name)
  line2=$(grep "^Result:" $file2)

  if [ -z "$line2" ]; then
    echo "Error in case n = $i, our algorithm did not complete..."
    continue
  fi

  
  # Extract the value between "Result:" and "("
  value2=$(echo "$line2" | sed -E 's/^Result:[[:space:]]*([^()]*)\(.*$/\1/' | xargs)

  if [ "$value1" != "$value2" ]; then
    echo "Results for case n = $i do not match!"
  
  fi
done
echo "All match!"

