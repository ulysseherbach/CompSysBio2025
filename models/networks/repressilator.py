"""Basic repressilator network (3 genes)."""
from models.networks import Network

# Initialize network
network = Network(3)

# Set basal parameters
network.basal[:] = 5

# Set interaction matrix
network.inter[0, 1] = -10
network.inter[1, 2] = -10
network.inter[2, 0] = -10
