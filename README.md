# Your Project

A short description of what your project does.  

## Installation

Install in editable mode with dev dependencies (so you get pytest, linting, sphinx docs):

```bash
pip install -e ".[dev]"


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