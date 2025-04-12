#!/bin/bash 

for i in $(seq 10 10);
do 

  { echo -e "Generating uba-$i \n"; } > "uba-$i.txt"
  python3 uba-gen.py $i

  #UBA model checking !use timeout on LINUX 
  gtimeout 30m ../../prism/bin/prism -explicit -ltluba -ubaverbosity 1 random2.pm -pf 'P=?[ HOA: {"'"UBAs/uba-$i.hoa"'", "sigma" <- s=0, "pi" <- s=1, "hash" <- s=2, "dollar" <- s=3}]' >> "uba-$i.txt"

done
