#### Build project:
```console
cd prism
make
```
#### Run simple examples
Model checking for UBA and LMC:

```console
bin/prism -explicit -ltluba -ubaverbosity 4 ../case-studies/family-of-UBAs/random2.pm -pf 'P=?[ HOA: {"UBAs/uba-3.hoa", "sigma" <- s=0, "pi" <- s=1, "hash" <- s=2, "dollar" <- s=3}]'
```

Model checking for GFG and LMC:

```console
bin/prism -explicit -gfgmc -ubaverbosity 4 ../case-studies/family-of-UBAs/random2.pm -pf 'P=?[ HOA: {"GFGs/gfg-3.hoa", "sigma_0" <- s=0, "pi_0" <- s=1, "hash_0" <- s=2, "dollar_0" <- s=3 }]'
```



## Run Experiments
#### Case study 1 - family of UBAs in the paper

Run experiments to compare the UBA model checking algorithm and the GFG model checking algorithm:
```console
cd ../case-studies
./run.sh
```

The LMC (in random2.pm) with 4 states generates the first label (sigma) with probablity one and uniform randomly generates the four labels afterwards.

cutoff time = 30 mins

We report the number of states of the automata, the number of states of the product and the running time.
|         | UBA   | GFG   | 
| :---:   | :---: | :---: |
| k=3 | #states = 7, #product = 20, time = 0.122 s | #states = 5, #product = 19, time = 0.205 s   |
| k=4 | #states = 15, #product = 40, time = 0.13 s   | #states = 6, #product = 22, time = 0.224 s   |
| k=5 | #states = 31, #product = 80, time = 0.228 s   | #states = 7, #product = 25, time = 0.222 s   |
| k=6 | #states = 63, #product = 160, time = 0.485 s   | #states = 8, #product = 28, time = 0.246 s   |
| k=7 | #states = 127, #product = 320, time = 2.873 s   | #states = 9, #product = 31, time = 0.247 s   |
| k=8 | #states = 255, #product = 640, time = 34.586 s   | #states = 10, #product = 34, time = 0.284 s   |
| k=9 | #states = 511, #product = 1280, time = 481.043 s| #states = 11, #product = 37, time = 0.309 s   |
| k=10 |#states = 1023, #product = 2560, timeout   | #states = 12, #product = 40, time = 0.365 s   |
