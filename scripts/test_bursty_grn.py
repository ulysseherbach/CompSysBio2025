"""Basic repressilator network (3 genes)."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from models.networks import Network
from models.solution import BurstyGRN, LimitGRN

# Example network
network = Network(3)
network.basal[:] = 5
network.inter[0, 1] = -10
network.inter[1, 2] = -10
network.inter[2, 0] = -10

model = BurstyGRN(network)
model.degradation_rate = 1
model.burst_frequency_max = 500
model.burst_size = model.degradation_rate/model.burst_frequency_max

# Time points
time = np.linspace(0, 20, 1000)

# Simulation of the PDMP model
sim = model.simulate(time, verb=True)

# Simulation of the ODE model (slow-fast limit)
model_ode = LimitGRN(network)
model_ode.degradation_rate = model.degradation_rate
model_ode.burst_frequency_max = model.burst_frequency_max
model_ode.burst_size = model.burst_size
sim_ode = model_ode.simulate(time, init_state=[0, 0.1, 0.2])

# Figure
fig = plt.figure(figsize=(10, 6))
gs = gridspec.GridSpec(2, 1, hspace=0.6)

# Plot protein levels
ax = plt.subplot(gs[0, 0])
for i in range(3):
    ax.set_title('Protein levels')
    ax.plot(sim.t, sim.x[:, i], label=f'$P_{{{i+1}}}$')
    ax.set_xlim(sim.t[0], sim.t[-1])
    ax.set_ylim(0, np.max([1.2*np.max(sim.x), 1]))
    ax.legend(loc='upper left', ncol=4, borderaxespad=0, frameon=False)

# Plot protein levels (ODE model)
ax = plt.subplot(gs[1, 0])
for i in range(3):
    ax.set_title(r'Protein levels - ODE model ($d_0/d_1\to\infty$)')
    ax.plot(sim_ode.t, sim_ode.x[:, i], label=f'$P_{{{i+1}}}$')
    ax.set_xlim(sim_ode.t[0], sim_ode.t[-1])
    ax.set_ylim(0, 1)
    ax.legend(loc='upper left', ncol=4, borderaxespad=0, frameon=False)

fig.savefig('bursty_grn.pdf', bbox_inches='tight')
