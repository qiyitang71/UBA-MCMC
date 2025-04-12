#!/bin/bash
#!/bin/bash

# Folder containing the files
results_folder="results"

# Loop through files gfg-3.txt to gfg-12.txt
for i in {3..12}; do
  file1="$results_folder/gfg-$i.txt"
  file2="$results_folder/uba-$i.txt"

  if [[ -f $fileteams1 ]]; then
    if [[ $i -gt 8 ]];then
         echo '\rowcolor{green!30!yellow!30}' 
    fi
    echo "$i &"
    # Add your processing commands below
    ./get_uba_latex.sh $file2
    echo "&"
    ./get_gfg_latex.sh $file1
    echo '\\\hline'
  else
    echo "File $file1 not found"
  fi
done

echo "Processing completed."
