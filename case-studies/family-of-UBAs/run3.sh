#!/bin/bash 

for i in $(seq 10 20);
do 
  { echo -e "Generating uba-$i ..."; } 2>&1 | tee "results/uba-$i.txt" 
  python3 uba-gen.py $i
  
  if [ "$i" -lt 20 ]; then
    #UBA model checking  
    echo -e "UBA Model Checking uba-$i ..."
    ../../prism/bin/prism -javamaxmem 100g -javastack 1g -timeout 30m -explicit -ltluba -ubaverbosity  1 -ubapower random3.pm -pf  'P=?[ HOA: {"'"UBAs/uba-$i.hoa"'", "sigma" <- "sigma_0", "pi" <- "pi_0", "hash" <- "hash_0", "dollar" <- "dollar_0"}]'  >> "results/uba-$i.txt"
  fi

  #timing the transformations
  { echo -e "Transforming to GFG-NCA gfg-$i ..."; } 2>&1 | tee "results/gfg-$i.txt"
  TIMEFORMAT=%R
  { time ../../../ltl2pba/uba2pba -a "UBAs/uba-$i.hoa" -o "GFGs/gfg-$i.hoa"; } >> "results/gfg-$i.txt" 2>&1

  { echo -e ""; } >> "results/gfg-$i.txt"
  #GFG model checking !use timeout on LINUX
  echo -e "GFG Model Checking gfg-$i ..." 
  ../../prism/bin/prism -javamaxmem 100g -javastack 1g -timeout 30m -explicit -gfgmc -ubaverbosity 1 -gfgpower random3.pm -pf 'P=?[ HOA: {"'"GFGs/gfg-$i.hoa"'", "sigma_0" <- "sigma_0", "pi_0" <- "pi_0", "hash_0" <- "hash_0","dollar_0" <- "dollar_0" }]' >> "results/gfg-$i.txt"


done
