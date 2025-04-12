#!/bin/bash 

#generating uba-3
echo -e "Generating UBA ..." 
python3 uba-gen.py 3

#timing the transformations
{ echo -e "Transforming and minimizing UBA ..."; }
TIMEFORMAT=%R
exec_time=$({ time ../uba2pba/uba2pba -a "UBAs/uba-3.hoa" -o "GFGs/gfg-3.hoa"; } 2>&1 > /dev/null)
echo "Transformation time = $exec_time"
