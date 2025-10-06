"""BurstyBase model (single gene with no feedback)."""
import numpy as np
import matplotlib.pyplot as plt
from models import BurstyBase
# from models.solution import BurstyBase

model = BurstyBase()
model.burst_size = 0.5
model.burst_frequency = 1
model.degradation_rate = 0.8

# Compute a single trajectory
time = np.linspace(0, 5, 1000)
sim = model.simulate(time, init_state=1)

# Show the trajectory
fig = plt.figure(figsize=(12, 3))
plt.plot(sim.t, sim.x)
plt.xlim(sim.t[0], sim.t[-1])
plt.ylim(0)
fig.savefig('bursty_base.pdf', bbox_inches='tight')
