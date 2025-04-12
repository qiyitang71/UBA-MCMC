#!/bin/bash

# Folder containing the files
results_folder="results"

# Format the row as a table with an extra "n" column
printf "| %-2s | %-10s | %-8s | %-8s | %-10s | %-8s | %-10s | %-10s | %-8s |\n" \
  "n" "UBA States" "Product" "Time" "GFG States" "Product" "Total Time" "Trans Time" "MC Time"

# Loop through files gfg-3.txt to gfg-12.txt and uba-3.txt to uba-12.txt
for i in {3..12}; do
  file1="$results_folder/gfg-$i.txt"
  file2="$results_folder/uba-$i.txt"

  # Run the scripts and capture the outputs
  uba_info=$(./get_uba_info.sh "$file2")
  gfg_info=$(./get_gfg_info.sh "$file1")

  # Extract the values using sed
  uba_states=$(echo "$uba_info" | sed -n 's/.*#states = \([0-9]*\).*/\1/p')
  uba_product=$(echo "$uba_info" | sed -n 's/.*#product = \([0-9]*\).*/\1/p')
  uba_time=$(echo "$uba_info" | sed -n 's/.*time = \([0-9]*\.[0-9]*\) s.*/\1/p')

  gfg_states=$(echo "$gfg_info" | sed -n 's/.*#states = \([0-9]*\).*/\1/p')
  gfg_product=$(echo "$gfg_info" | sed -n 's/.*#product = \([0-9]*\).*/\1/p')
  gfg_total_time=$(echo "$gfg_info" | sed -n 's/.*total_time = \([0-9]*\.[0-9]*\) s.*/\1/p')
  gfg_trans_time=$(echo "$gfg_info" | sed -n 's/.*: \([0-9]*\.[0-9]*\) (trans).*/\1/p')
  gfg_mc_time=$(echo "$gfg_info" | sed -n 's/.*, \([0-9]*\.[0-9]*\) (mc).*/\1/p')

  # Print the table row with the value of n (i)
  printf "| %-2s | %-10s | %-8s | %-8s | %-10s | %-8s | %-10s | %-10s | %-8s |\n" \
    "$i" "$uba_states" "$uba_product" "$uba_time" "$gfg_states" "$gfg_product" "$gfg_total_time" "$gfg_trans_time" "$gfg_mc_time"
done

echo "Table printed."
