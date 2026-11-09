# Your Project

Identifying dynamically influential nodes in biological networks is a central problem in systems biology, particularly for prioritizing intervention targets in gene regulatory networks. In this paper, we propose a \textcolor{blue}{target-specific} Shapley-value-based framework for assessing the importance of nodes in a Boolean network. The framework comprises two complementary measures: the \emph{Knock-out} and the \emph{Knock-in Shapley values}. \textcolor{blue}{Moreover, by exploiting the logical structure of the network, we present a propagation-based method to avoid exhaustive simulations}. The approach is exact for acyclic networks and provides good approximations for cyclic networks. Evaluation on benchmark models from the Cell Collective database shows that the propagation method accurately recovers node importance rankings while achieving substantial speed-ups. 

## Installation

Install in editable mode with dev dependencies (so you get pytest, linting, sphinx docs):

```bash
pip install -e ".[dev]"
```
## USAGE 
```bash
usage: __main__.py [-h] -n NAME -e EXPRESSION -o OUTPUT [-r SAMPLING] [-p SAMPLEPER] [-c SIMPER] [-k SCHEME] [-m MODE] [-s SIMPROP] [-v] [-d] [-t]

Input the expression file and

options:
  -h, --help            show this help message and exit
  -n NAME, --name NAME  Name of the system
  -e EXPRESSION, --expression EXPRESSION
                        path to expression file
  -o OUTPUT, --output OUTPUT
                        list of interested components to examine
  -r SAMPLING, --sampling SAMPLING
                        Sampling strategy for input states: "all" (default), "uniform", "antithetic", or "stratified"
  -p SAMPLEPER, --samplePer SAMPLEPER
                        Percentage of combinations to sample if using a sampling strategy different from "all". Default is 60.0
  -c SIMPER, --simPer SIMPER
                        Percentage of nodes to be simulated. Default is 0.0
  -k SCHEME, --scheme SCHEME
                        Scheme of calculation, can be "Shapley" or "Uniform", default is "Shapley"
  -m MODE, --mode MODE  Mode of calculation, can be "Test" or "Propagation" or "Simulation", default is "Propagation"
  -s SIMPROP, --simprop SIMPROP
                        Probability to simulate an original node during the propagation
  -v, --verbose         Print log or not
  -d, --decimal         Log decimal mappings for hardcore debug mode
  -t, --test            Perform full network perturbation to test with the obtained result
```
## BATCH TEST 
To run batch test with the benchmark models
python run_batch.py

### To choose the testing scheme 
Modify this part of the file run_batch.py:
```bash
# ---------------- CONFIG (hard-coded) ----------------
INPUT_LIST = "realbenchmark.txt"            # TAB-separated: path<TAB>out1<TAB>out2<TAB>out3
LOG_DIR    = Path("logs/realbenchmark_Prop_sample_60_5thRun")  # where the terminal logs are saved 
PROGRAM    = [sys.executable, "-m", "Shapley"]
ECHO_TO_TERMINAL = True                  # also print the log live

# Fixed parameters (edit freely; remove any you don't want)
FIXED_ARGS = [
    "-r", "uniform", # Sampling strategy: all, uniform, antithetic, stratified
    "-p", "60.0", # Sample percentage (0-100)
    "-c", "0.0", # Percentage of nodes to be simulated (0-100)
    "-k", "Shapley", # Weighting scheme: Shapley, Uniform
    "-m", "Propagation", # Method: Propagation, Simulation
    "-t" # Enable Test mode to compare with fully simulation results
]
```
with the args as from the usege of the main function.