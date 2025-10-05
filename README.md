# CompSysBio2025 hands-on session

*Simulation and reverse-engineering of mechanistic GRN-driven models of gene expression*

Material for hands-on session at [CompSysBio2025](https://project.inria.fr/compsysbio2025/) (Advanced Lecture Course on Computational Systems Biology).

This repository is dedicated to the first part. Here we deal with simulation of mechanistic models, including stochastic gene expression.

## Quickstart

1. Clone the repository (or just download the zip archive) and put the folder where you want
2. Optionally, open a terminal at this location and perform editable installation:  
    `pip install -e .`
3. Check if everything works by running the `main.py` script

Now you can start! The files to be modified are:

- `models/bursty_base.py`
- `models/bursty_grn.py`
- `models/ode_grn.py`

**Note:** If step 2 is not done, you will not be able to run directly the `models/...` files (which can be useful for testing). But you can always run the `main.py` script and load `models` from there.
