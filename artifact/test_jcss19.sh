#!/bin/bash 

#generating uba-3
echo -e "Generating UBA ..." 
python3 uba-gen.py 3

#UBA model checking using artifact JCSS19
echo -e "UBA Model Checking uba-3 ..." 
../prism/bin/prism -javamaxmem 200g -javastack 1g -timeout 30m -explicit -ltluba -ubaverbosity  1 -ubapower random_lmc.pm -pf  'P=?[ HOA: {"'"UBAs/uba-3.hoa"'", "sigma" <- "sigma", "pi" <- "pi", "hash" <- "hash", "dollar" <- "dollar"}]'  > "results/uba-3.txt"

#report results
echo -n "Output: "; ./get_uba_info.sh results/uba-3.txt
