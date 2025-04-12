#!/bin/bash 

for i in $(seq 3 10);
do 
  { echo -e "Generating dba-$i ..."; } 2>&1 | tee "results/dba-$i.txt" 
  python3 dba-gen.py $i

  #UBA model checking !use timeout on LINUX 
  echo -e "UBA Model Checking dba-$i ..."
  gtimeout 30m ../../prism/bin/prism -explicit -ltluba -ubaverbosity 1 random.pm -pf 'P=?[ HOA: {"'"DBAs/dba-$i.hoa"'", "sigma" <- s=0, "pi" <- s=1, "hash" <- s=2}]' >> "results/dba-$i.txt"

  #timing the transformations
  { echo -e "Transforming to GFG-NCA gfg-$i ..."; } 2>&1 | tee "results/gfg-$i.txt"
  TIMEFORMAT=%R
  { time ../../../ltl2pba/uba2pba -a "DBAs/dba-$i.hoa" -o "GFGs/gfg-$i.hoa"; } >> "results/gfg-$i.txt" 2>&1

  { echo -e ""; } >> "results/gfg-$i.txt"
  #GFG model checking !use timeout on LINUX
  echo -e "GFG Model Checking gfg-$i ..." 
  gtimeout 30m ../../prism/bin/prism -explicit -gfgmc -ubaverbosity 1 random.pm -pf 'P=?[ HOA: {"'"GFGs/gfg-$i.hoa"'", "sigma" <- s=0, "pi" <- s=1, "hash" <- s=2 }]' >> "results/gfg-$i.txt"

done
