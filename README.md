# _evosim_: a Python package for fast and scalable stochastic simulations of evolutionary dynamics

_evosim_ is a Python implementation of a fast, scalable, and unbiased tau-leaping algorithm for the simulation of continuous-time Markov birth-death-mutation branching processes. These stochastic processes are commonly simulated by Gillespie's Stochastic Simulation Algorithm (SSA), but the SSA does not scale up well to large population times and becomes computationally prohibitive at populations above a million cells. Many biomedical scenarios (e.g. simulation of tumor evolution and progression) can evolve cell populations that are much larger than this, calling for faster methods that retain simulation accuracy.

_evosim_ is a wrapper for the tau-leaping algorithm, described in detail in the Supplemental Information in **PAPER**. The package provides flexibility in terms of changing evolutionary parameters over time (e.g. due to changes in evolutionary conditions, like changes in growth due to drug administration or changes in mutation rates). The package provides tracking and plotting capabilities, a function for subsampling from the entire population (e.g. as in an serial dilution or surgical resection), as well as functions for computing a number of diversity indices (e.g. Shannon and Gini indices) and for converting the time units on the output.

## Getting Started

### Dependencies

* Python>=3.8
* numpy>=1.20.0
* pandas>=1.2.4
* matplotlib>=3.3.4
* seaborn>=0.11.1


### Installing

To install from GitHub: 

git clone https://github.com/daliten/evosim

pip install ./evosim

## License

_evosim_ is released under the General Public License version 3 (GPLv3+).

## Using evosim

_evosim_ contains two modules, `Simulator` for single-run simulations and `EnsembleSimulator` for multi-run simulations. They are launched as
```
from evosim import Simulator, EnsembleSimulator

sim=Simulator()
esim=EnsembleSimulator()
```
and each come with appropriate functionalities that include tracking and plotting the clonal distribution over time and computing diversity indices (single runs) and distributions of diversity indices (ensemble). Both simulators are constructed to allow a series of simulation steps to be performed sequentially so that evolutionary conditions (e.g. presence/concentration of drugs) can be changed with full flexibility throughout the course of the simulation. See Directory for examples and tutorials and Documentation.

Both types of simulators can be run either in "leaping" mode (default), employing the fast and scalable unbiased tau-leaping algorithm presented in **PAPER**, or in "SSA" mode, which uses the SSA without any modifications.

The parameters that must be supplied to either simulator are the clonal distribution at the simulation start (species/phenotypes names or labels with respective cell numbers), mutation probabilities i->j (probability of mutation to phenotype j in a single division of a phenotype i cell), birth/division rates (inverse expected time, in user-chosen time units, to the division of an individual cell of each phenotype), death rates (inverse expected time, in user-chosen time units, to the death of an individual cell of each phenotype), overall simulation time (in user-chosen time units), and the simulation interval for recording sequential population distributions (in user-chosen time units). Please note that any time units may be used in the simulation; however, this choice must be consistent across all parameters. _evosim_ plots provide an option for time conversion from/to different units (e.g. if user specified parameters are provided in days, the time course of the clonal distribution can be plotted in weeks, months, or years instead).

## Documentation

See Documentation for full functional documentation and **PAPER** for background, algorithm, and pseudocode.

## Examples and tutorials

See **examples_tutorials** for examples and tutorials. The EnsembleSimulator example, which generates, treats, and analyzes a cancer patient cohort, assumes familiarity with the functionalities described in Simulator example 1. 
