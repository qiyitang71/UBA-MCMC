#!/bin/bash 

#generating gfg-3
echo -e "Generating GfG automaton ..." 
python3 gfg-gen.py 3

#GFG model checking
echo -e "GFG Model Checking gfg-3 ..." 
../prism/bin/prism -javamaxmem 100g -javastack 1g -timeout 30m -explicit -gfgmc -ubaverbosity 1 -gfgpower random_lmc.pm -pf 'P=?[ HOA: {"'"GFGs/gfg-3.hoa"'", "sigma_0" <- "sigma", "pi_0" <- "pi", "hash_0" <- "hash","dollar_0" <- "dollar" }]' > "results/gfg-3.txt"

#report results
echo -n "Output: "; ./get_gfg_info.sh results/gfg-3.txt

