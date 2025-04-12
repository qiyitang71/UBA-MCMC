# Artifact for Paper 151

This artifact provides the source code for our tool, [UBA-MCMC](https://github.com/qiyitang71/UBA-MCMC), which performs model checking of Markov chains against unambiguous Büchi automata (UBA) specifications.  
UBA-MCMC is distributed under the GNU General Public License (GPL) v2.

## Load Docker Container

Ensure [Docker](https://www.docker.com/get-started/) is installed.  
**Note:** Docker commands typically require **root** access.

Load the Docker image:

```bash
sudo docker load < uba-mcmc.tar.gz
# or
sudo docker load -i uba-mcmc.tar.gz
```

Check the image has been loaded:

```bash
sudo docker images
```

Run the container:

```bash
sudo docker run -it uba-mcmc
```

## Benchmarks

The random labelled Markov chain used to reproduce the results is:

```
./artifact/random_lmc.pm
```

The set of benchmarks (UBAs) are generated on the fly and stored in:

```
./artifact/UBAs/
```

### File Input Format

All automata are written in the [HOA format](https://adl.github.io/hoaf/).  
An example automaton is shown below:

```hoa
HOA: v1
name: "FGa"
States: 4
Start: 0
AP: 1 "a"
acc-name: Buchi
Acceptance: 1 Inf(0)
properties: trans-labels explicit-labels state-acc complete unambiguous
properties: stutter-invariant
--BODY--
State: 0
[!0] 0
[0] 1
[0] 2
State: 1 {0}
[0] 1
[!0] 3
State: 2
[!0] 0
[0] 2
State: 3
[t] 3
--END--
```

## Installation

Our tool includes two components: automata transformation (`uba2pba`) and GfG model checking (`prism-gfg`).

First, compile `uba2pba`:

```bash
cd uba2pba
make
```

Then compile `prism-gfg`:

```bash
cd ../prism
make
```

## Smoke Tests

We compare our method with the state-of-the-art JCSS19 artifact ([link](https://wwwtcs.inf.tu-dresden.de/ALGI/TR/JCSS19/)).

Go to the `artifact` directory:

```bash
cd ../artifact
```

1. Test JCSS19:

```bash
./test_jcss19.sh
```

Expected output (time may vary):

```
Generating UBA ...
UBA Model Checking uba-3 ...
Output: #states = 7, #product = 3779, time = 0.327 s
```

2. Test `uba2pba`:

```bash
./test_transformation.sh
```

This produces `GFGs/gfg-3.hoa`.

3. Test GfG model checking:

```bash
./test_gfg.sh
```

Expected output (time may vary):

```
Generating GfG automaton ...
GFG Model Checking gfg-3 ...
Output: #states = 5, #product = 3910, total_time = 1.182 s: 0 (trans), 1.182 (mc)
```

## Full Experiments

To reproduce the full benchmark results:

```bash
./run.sh
```

This may take around **2 hours**. Results are saved in `./results`.

To print the results:

```bash
./print_data.sh
```

The results can be compared with Table 1 in the paper. While timings may vary, the trend remains: from `n = 8`, our algorithm outperforms JCSS19.  
Due to memory demands, JCSS19 may run out of resources before our method encounters issues.

To exit Docker:

```bash
exit
```

## Model Checking Another Random Markov Chain

Generate a new Markov chain:

```bash
python3 dtmc-gen.py lmc-test.pm
```

Run model checking:

```bash
../prism/bin/prism -javamaxmem 100g -javastack 1g -timeout 30m -explicit -gfgmc -ubaverbosity 1 -gfgpower lmc-test.pm -pf 'P=?[ HOA: {"'"GFGs/gfg-3.hoa"'", "sigma_0" <- "sigma", "pi_0" <- "pi", "hash_0" <- "hash","dollar_0" <- "dollar" }]'
```

## Build Docker Container

Clone the repo:

```bash
git clone https://github.com/qiyitang71/UBA-MCMC.git
```

Build the image:

```bash
cd UBA-MCMC
sudo docker build -t uba-mcmc .
```

## Remove Docker Container

List containers:

```bash
sudo docker ps -a
```

Remove a container:

```bash
sudo docker rm <container_id>
```
