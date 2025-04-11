# uba2pba

`uba2pba` is a tool developed as part of the ongoing ltl2pba project. It takes an unambiguous Büchi automaton (UBA) as input, determinises it by making nondeterministic choices explicit, and then minimises the resulting deterministic automaton using the techniques described in [1].



### Requirements
* [Spot model checker](https://spot.lrde.epita.fr/)

To install Spot manually:
```
./configure
make && make install
```

**Note**: in our provided virtual machine, the Spot library is already pre-installed.

### Compilation
To compile `uba2pba`, simply run:
```
make
```

This will generate an executable file named **uba2pba** !

### Usage

Run `uba2pba` with an input UBA stored in "filename": 
```
./uba2pba -a filename
```
By default, the output automaton is printed to standard output!

To save the result automaton to a file, use:
```
./uba2pba -a filename -o out_filename
```

To check whether the result is language-equivalent to the determinised automaton, use:
```
./uba2pba -a filename -o out_filename -c
```
A successful equivalence check will print
```
Equivalence: 1
```
indicating that the output automaton is correct.

### Smoke test
As a quick test, consider a deterministic Büchi automaton (DBA) in `dba.hoa` with 6 states. 
Run:
```
./uba2pba -a dba.hoa
```
This should produce a minimised automaton with 4 states.

To save the resultant automaton, e.g. to `gfg.hoa` and verify its correctness , run:
```
./uba2pba -a dba.hoa -o gfg.hoa -c
```
You are expected to see the following:
```
The UBA is deterministic
The max outgoing degree is 0
Minimising DCA
UBA has 7 states
Computing safe components
Computing language equivalence
Now normalise DCA
Computing subsafe language equivalence
Now safe centralise DCA
Now safe minimise DCA
Now output DCAs
GfG-NCA has 4 states
Equivalence: 1
```
The final line confirms that the resultant automaton is correct.




---
[1] Bader Abu Radi, Orna Kupferman:
Minimization and Canonization of GFG Transition-Based Automata. Log. Methods Comput. Sci. 18(3) (2022)