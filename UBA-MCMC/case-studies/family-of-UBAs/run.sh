#!/bin/bash 

for i in $(seq 3 10);
do 
  python3 uba-gen.py $i
  gtimeout 30m ../../prism/bin/prism -explicit -ltluba -ubaverbosity 1 random2.pm -pf 'P=?[ HOA: {"'"UBAs/uba-$i.hoa"'", "sigma" <- s=0, "pi" <- s=1, "hash" <- s=2, "dollar" <- s=3}]' > results/"uba-$i.txt"
  #use timeout on LINUX 
  python3 gfg-gen.py $i
  gtimeout 30m ../../prism/bin/prism -explicit -gfgmc -ubaverbosity 1 random2.pm -pf 'P=?[ HOA: {"'"GFGs/gfg-$i.hoa"'", "sigma_0" <- s=0, "pi_0" <- s=1, "hash_0" <- s=2, "dollar_0" <- s=3 }]' > results/"gfg-$i.txt"
done

../../prism/bin/prism -explicit -gfgmc -ubaverbosity 4 random2.pm -pf 'P=?[ HOA: {"'"GFGs/gfg-3.hoa"'", "sigma_0" <- s=0, "pi_0" <- s=1, "hash_0" <- s=2, "dollar_0" <- s=3 }]'