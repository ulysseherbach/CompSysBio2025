"""CompSysBio2025 hands-on session.

Here we simulate mechanistic models of gene expression driven by gene
regulatory networks (GRNs), with an emphasis on piecewise-deterministic
Markov processes (PDMPs) to describe biological stochasticity at the
single-cell level.
"""
from importlib.metadata import version as _version
from models.bursty_base import BurstyBase
from models.bursty_grn import BurstyGRN
from models.limit_grn import LimitGRN

__all__ = ['BurstyBase', 'BurstyGRN', 'LimitGRN']

try:
    __version__ = _version('models')
except Exception:
    __version__ = 'unknown version'
