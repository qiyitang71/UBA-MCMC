
Build project:
```console
cd prism
make
```

Model checking for UBA and LMC:

```console
bin/prism -explicit -ltluba -ubaverbosity 4 ../prism-examples/dice/random2.pm -pf 'P=?[ HOA: {"uba-family-state-acc-3.hoa", "sigma" <- s=0, "pi" <- s=1, "hash" <- s=2, "dollar" <- s=3}]'
```

Model checking for GFG and LMC:

```console
bin/prism -explicit -gfgmc -ubaverbosity 4 ../prism-examples/dice/random2.pm -pf 'P=?[ HOA: {"gfg-family-trans-acc-3.hoa", "sigma_0" <- s=0, "pi_0" <- s=1, "hash_0" <- s=2, "dollar_0" <- s=3 }]'
```

