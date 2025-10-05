# CompSysBio2025 hands-on session

*Simulation and reverse-engineering of mechanistic GRN-driven models of gene expression*

This repository is dedicated to the first part: here we deal with simulating mechanistic models of gene expression driven by gene regulatory networks (GRNs), with an emphasis on piecewise-deterministic Markov processes (PDMPs) to describe biological stochasticity at the single-cell level.

## Quickstart

1. Clone the repository (or just download the zip archive) and put the folder where you want
2. Optionally, open a terminal at this location and perform editable installation:  
    `pip install -e .`
3. Check if everything works by running the `main.py` script

Now you can start! The files to be modified are:

- `models/bursty_base.py`
- `models/bursty_grn.py`
- `models/limit_grn.py`

**Note:** If step 2 is not done, you will not be able to run directly the `models/...` files (which can be useful for testing). But you can always run the `main.py` script and load `models` from there.
