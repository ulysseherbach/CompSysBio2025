"""Useful functions for defining interactions."""
import numpy as np
from scipy.special import expit


class Network:
    """Store network parameters (`basal` and `inter`)."""

    def __init__(self, n_genes):
        self.n_genes = n_genes  # Number of genes
        self.basal = np.zeros(n_genes)  # Basal activities
        self.inter = np.zeros((n_genes, n_genes))  # Interaction matrix


def kon_sigmoid(x, k0, k1, basal, inter):
    """Define interactions using a logistic function (sigmoid)."""
    sigma = expit(basal + x @ inter)
    return (1-sigma)*k0 + sigma*k1
