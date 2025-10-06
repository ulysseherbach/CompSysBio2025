"""Basic toggle switch network (2 genes)."""
from models.networks import Network

# Initialize network
network = Network(2)

# Set basal parameters
network.basal[:] = 5

# Set interaction matrix
network.inter[0, 1] = -10
network.inter[1, 0] = -10
