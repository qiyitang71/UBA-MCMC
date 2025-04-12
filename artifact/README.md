# Artifact for Paper 108

In this artifact, you can find the souce code of our tool [DFAMiner](https://github.com/liyong31/DFAMiner.git) that learns the minimal separating DFAs from labelled samples.
DFAMiner is distributed under the MIT license.

## User info
Both the user name and the password of the Ubuntu system are `test`

## Benchmarks
The set of benchmarks we used for the experiments are in [dataset](./dataset/) directory.
Each index/directory N in [dataset](./dataset/) has 100 benchmark instances of random samples, where in each benchmark instance, there are 50 × N samples.
The alphabet for the samples has two symbols.

#### File input format

All tools accept a file in [Abbadingo format](https://abbadingo.cs.nuim.ie/).
An example file is given below:
```
8 2
1 3 0 0 0
1 3 0 0 1
0 3 0 1 0
0 3 0 1 1
1 3 1 0 0
0 3 1 0 1
0 3 1 1 0
0 3 1 1 1
```
The first line gives the number of samples, i.e., 8 and the size of the alphabet, i.e., 2.
An alphabet of size n will have letters from {0, ..., n-1}.
Each line after the first line will first specify the membership of the word, then the length of the word and afterwards the word sequence.
Here `1` means accept, `0` reject and `-1` don't care.
We also call words labelled with `1` the *positive* words, and words labelled with `0` the *negative* words.
For instance, the second line `1 3 0 0 0` gives a positive word `000` of length 3.


## Tools
We compare with two state of the art tools, [DFAInductor](https://github.com/ctlab/DFA-Inductor-py) and [DFAIdentify](https://github.com/mvcisback/dfa-identify).


#### Tests

To test whether all tools have been installed correctly, we can do the following under the **artifact** directory.

1. Test DFAMiner:
```
python3 ./DFAMiner/dfaminer.py --file example.txt --out example.dot 
```
DFAMiner takes as input the sample file *example.txt* and outputs the learned minimal separating DFA to *example.dot* in [DOT](https://graphviz.org/doc/info/lang.html) format.
We can use xdot to visualise the resultant automaton.
```
xdot example.dot
```  

2. Test DFAInductor:
```
dfainductor -i example.txt -o example-ind.dot -s cadical153
```
For DFAInductor, we also need to specify which SAT solver you wish DFAInductor to use.
By default, DFAMiner and DFAIdentify will use CaDiCal 1.5.3.
We refer to [PySAT](https://pysathq.github.io/features/) homepage for more available SAT solvers.  

3. Test DFAIdentify:
```
python3 ./dfa-identify/dfaid.py --file example.txt --out example-id.dot
```
One can verify that all learned DFAs accept all positive words and reject all negative words.
The learned DFAs by all tools should be the same, since for our example, there is only one possible minimal separating DFA.
In general, the learned minimal DFAs may be different, but they all should have the same number of states.

## Experiments

### Original data
In the `results/orginal` folder, you can find all log files of the experiments we conducted for DFAMiner and DFAInductor in the paper.
To replicate part of the Table 3, go to directory `results/orginal` and type the following commands
```
python3 gen_table_3.py
```
One can see the output and the corresponding data in Table 3; they should be the same.
Note that we only included original data for DFAInductor and DFAMiner.

### Replication of experiments

To run all the experiments reported in this paper, use the following command:
```
```
However, this may take more than **20 days** to complete the experiments on an Intel i7-483 4790 3.60 GHz processor with 16GB of RAM.

Instead, one can choose to run the experiments up to index N=10 using the following command.
```
```
This set of benchmarks should take around **6 hours** to complete.


## Applications
In this application, we try to extract the minimal reachability or safety automata for solving parity games.
The idea is to extract a minimal reachability/safety automaton that accepts all *even* cycles and reject all *odd* cycles.
A cycle is a word fragment of length greater than 1 in which the first letter is the same as the last letter.
A cycle is even (respectively, odd), if the maximal number occurs inside the cycle is even (respectively, odd).

For instance, the word $(001212)^{\omega}$ contains only even cycles 00,121 and 212, while the word $(1312331)^{\omega}$ contains only odd cycles such as 131, 3123, and 33.
It is shown in [1,2] that given a parity game G, and a deterministic reachability/safety automaton D that accepts only even cycles and reject odd cycles, we can reduce the parity game to the reachability/safety game G×D.
For more details, we refer to [1,2].

Although the length of words in the samples should be infinite, but we can also use finite words of a given length, as we are looking for reachability/safety DFAs.
As long as the length is large enough, the DFA will converge to the minimal one over the infinite words.

We apply our tool DFAMiner to extract the minimal reachability/safety DFA from the set of *positive* words that only contain even cycles and the set of *negative* words that contain only odd cycles.
For a given number of colours, say n, and the length of the word m, the number of words that need to generated is $m^n$.
For instance, when the number of colours is 2 and the length is 3, all the possible samples can be found in our example.txt file.

If the reviewer wish to generate the sample file for colour 4 and length 11 (see Table 2), use the following command:
```
python3 scripts/paritygen.py --c 4 --l 7 --out parity/data4-7.txt
```
The output should be the following:
```
Completed generation of samples...
#colour=4
#length=7
#pos_words=1645
#neg_words=5235
#dont_words=9504
```
This means we have generated 1645 positive words, 5235 negative words and 9604 don't-care words.
The don't-care words contain both even and odd cycles.

To learn the minimal DFA with DFAMiner, use the following command:
```
python3 DFAMiner/dfaminer.py --file parity/data4-7.txt --out min-dfa4-7.dot --safety
``` 
It will take less than one second to output the minimal DFA to file min-dfa4-7.dot.

DFAInductor will take around 3 seconds to give the output.
```
dfainductor -i parity/data4-7.txt -o dfaind-min-dfa4-7.dot -s cadical153
```

When n is 6, the length of the words needs to reach 15, in order to allow the DFA converge to the target.
That is, we need to generate $6^{15}$ words, which is far too many for the artifact. (The sample file would take up more than 10GB of space.)
Therefore, we do not recommend the reviewer to generate the sample file for the number of colours up to 4.
Nonetheless, we still provide the constructed 3DFA for number of colours greater than 4 in `./parity` folder.
In order to allow DFAMiner to take as input the 3DFAs, use the tool minimser in DFAMiner instead:
```
python3 DFAMiner/minimiser.py --f parity/dfa5-11.txt --out min-dfa5-11.dot --safety
```

To replicate Table 2, first generate all samples:
```
cd scripts
chmod +x gen_samples.sh
./gen_samples.sh
cd ..
```
One can see that the number of samples should be consistent with those in Table 2 for colours 2, 3, and 4.
To reproduce the numer of states in the minimal DFA, run the following commands and check for each case:
```
python3 DFAMiner/dfaminer.py --file parity/data2-3.txt --out parity/min-dfa2-3.dot --safety
python3 DFAMiner/dfaminer.py --file parity/data3-5.txt --out parity/min-dfa3-5.dot --safety
python3 DFAMiner/dfaminer.py --file parity/data4-7.txt --out parity/min-dfa4-7.dot --safety
python3 DFAMiner/minimiser.py --f parity/dfa5-11.txt --out min-dfa5-11.dot --safety
python3 DFAMiner/minimiser.py --f parity/dfa6-15.txt --out min-dfa6-15.dot --safety
```
One can also check whether increasing the length of the word help to get smaller DFA by the following command for colour 5:
```
python3 DFAMiner/minimiser.py --f parity/dfa5-12.txt --out min-dfa5-12.dot --safety
```
The minimal DFA still has 5 states, as before. This means we have reached a fixed point.


[1] Bojańczyk, M., Czerwiński, W.: An Automata Toolbox.  (2018)

[2] Czerwinski, W., Daviaud, L., Fijalkow, N., Jurdzinski, M., Lazic, R., Parys, P.: Universal Trees Grow Inside Separating Automata: Quasi-Polynomial Lower Bounds for Parity Games. In: Symposium on Discrete Algorithms 19. p. 2333–2349. SIAM (2019)




