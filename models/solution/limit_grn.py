"""Simulation of the limit model for a gene regulatory network."""
import numpy as np
import models._utils as utils
from models.networks import Network, kon_sigmoid


class Trajectory:
    """Store single-cell gene expression trajectories."""

    def __init__(self, t, x):
        self.t = t  # Time points
        self.x = x  # Protein levels


class LimitGRN:
    """Limit model for a GRN with any number of interacting genes."""

    def __init__(self, network: Network,
        burst_size=1.0,
        burst_frequency_min=0.0,
        burst_frequency_max=2.0,
        degradation_rate=1.0):

        # Set model parameters
        self.network = network
        self.burst_size = burst_size
        self.burst_frequency_min = burst_frequency_min
        self.burst_frequency_max = burst_frequency_max
        self.degradation_rate = degradation_rate

        # Store the number of genes
        self.n_genes = network.basal.size

    def kon(self, x):
        """Define burst frequencies as a function of protein levels."""
        k0 = self.burst_frequency_min
        k1 = self.burst_frequency_max
        basal = self.network.basal
        inter = self.network.inter
        return kon_sigmoid(x, k0, k1, basal, inter)

    def euler_step(self, dt, x):
        """Perform basic Euler step for the limit model."""
        burst_size = self.burst_size
        degradation_rate = self.degradation_rate
        return (1 - dt*degradation_rate)*x + dt*burst_size*self.kon(x)

    def simulate(self, time, init_state=None, verb=False):
        """Perform basic simulation (extracted at given time points)."""
        if init_state is None:
            init_state = np.zeros(self.n_genes)

        # Check simulation parameters
        time = utils.check_time_points(time).reshape((-1,))
        init_state = utils.check_init_state(init_state, shape=(self.n_genes,))

        # Set Euler step size
        dt = 1e-3 / self.degradation_rate

        # Initialize trajectory array
        traj = np.zeros((time.size, self.n_genes))

        # Initialize current time and state
        t, x = 0, init_state

        # Optional: record number of steps
        c = 0

        # Core loop for simulation and recording
        for k in range(time.size):
            while t < time[k]:

                x = self.euler_step(dt, x)
                t += dt

                # Optional: record number of steps
                c += 1

            # Record protein levels
            traj[k] = x

        # Display info about steps
        if verb:
            print(f'ODE simulation used {c} steps (step size = {dt:.5f})')

        return Trajectory(time, traj)


# Tests
if __name__ == '__main__':
    from models.networks import toggle_switch
    model = LimitGRN(toggle_switch)
    # Simulation
    time = np.linspace(0, 10, 3)
    init_state = [1, 0]
    sim = model.simulate(time, init_state, verb=True)
    print(sim.t)
    print(sim.x)
