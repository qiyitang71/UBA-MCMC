#!/bin/bash 

for i in $(seq 3 10);
do 

  { echo -e "Generating uba-$i \n"; } > results/"uba-$i.txt"
  python3 uba-gen.py $i
  #UBA model checking !use timeout on LINUX 
  gtimeout 30m ../../prism/bin/prism -explicit -ltluba -ubaverbosity 1 random2.pm -pf 'P=?[ HOA: {"'"UBAs/uba-$i.hoa"'", "sigma" <- s=0, "pi" <- s=1, "hash" <- s=2, "dollar" <- s=3}]' >> results/"uba-$i.txt"

  #timing the transformation
  { echo -e "Transforming to GFG-NCA gfg-$i \n"; } > results/"gfg-$i.txt"
  { time ../../../ltl2pba/uba2pba -a "UBAs/uba-$i.hoa" -o "GFGs/gfg-$i.hoa"; } >> results/"gfg-$i.txt" 2>&1
  { echo -e "\n"; } >> results/"gfg-$i.txt"
  #GFG model checking !use timeout on LINUX 
  gtimeout 30m ../../prism/bin/prism -explicit -gfgmc -ubaverbosity 1 random2.pm -pf 'P=?[ HOA: {"'"GFGs/gfg-$i.hoa"'", "sigma_0" <- s=0, "pi_0" <- s=1, "hash_0" <- s=2, "dollar_0" <- s=3 }]' >> results/"gfg-$i.txt"

done
