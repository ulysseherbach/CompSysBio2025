"""Main file for running tests."""
import numpy as np
from models.networks import repressilator
import models
# # Uncomment to use the solution
# import models.solution as models

# # 1. Single gene with no feedback
# model = models.BurstyBase()

# 2. Gene regulatory network
model = models.BurstyGRN(repressilator)

# # 3. Deterministic limit model
# model = models.LimitGRN(repressilator)

# Define time points to record
time = np.linspace(0, 50, 6)

# Simulate the model
sim = model.simulate(time, verb=True)

print(sim.t)
print(sim.x)
