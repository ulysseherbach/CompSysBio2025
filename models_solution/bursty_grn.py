"""Simulation of the Bursty model for a gene regulatory network."""
import numpy as np
import models._utils as utils
from models.networks import Network, kon_sigmoid


class Trajectory:
    """Store single-cell gene expression trajectories."""

    def __init__(self, t, x):
        self.t = t  # Time points
        self.x = x  # Protein levels


class BurstyGRN:
    """Bursty model for a GRN with any number of interacting genes."""

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

    def rate_bound(self):
        """Compute the global burst rate upper bound."""
        return self.n_genes * self.burst_frequency_max

    def flow(self, time, x):
        """Define the deterministic flow between jumps."""
        degradation_rate = self.degradation_rate
        # Explicit solution of the ODE part
        return x * np.exp(- degradation_rate * time)

    def random_step(self, x, rng: np.random.Generator):
        """Compute next jump waiting time and state just after jump."""
        tau = self.rate_bound()

        # Sample waiting time before next jump
        u = rng.exponential(scale=1/tau)

        # Update state just before jump
        x = self.flow(u, x)

        # Construct the vector of jump probabilities
        v = np.zeros(self.n_genes + 1)
        v[1:] = self.kon(x)/tau  # i = 1, 2, ... : burst of protein i
        v[0] = 1 - np.sum(v[1:])  # i = 0 : no change (phantom jump)

        # Sample from this distribution
        # NB: This implementation deals robustly with precision errors
        i = np.searchsorted(np.cumsum(v), rng.random(), side='right')

        # Test if jump is real (i > 0) or phantom (i = 0)
        is_jump = i > 0

        # Perform the jump
        if is_jump:
            x[i-1] += rng.exponential(self.burst_size)

        return u, x, is_jump

    def simulate(self, time, init_state=None, seed=None, verb=False):
        """Perform exact simulation (extracted at given time points)."""
        if init_state is None:
            init_state = np.zeros(self.n_genes)

        # Check simulation parameters
        time = utils.check_time_points(time).reshape((-1,))
        init_state = utils.check_init_state(init_state, shape=(self.n_genes,))

        # Define a random generator
        rng = np.random.default_rng(seed)

        # Optional: record jump counts (phantom and true)
        n_jumps = np.zeros(2, dtype=np.uint)

        # Initialize trajectory array
        traj = np.zeros((time.size, self.n_genes))

        # Initialize current time and state
        t, x = 0, init_state

        # Initialize previous time and state
        t_old, x_old = t, x

        # Core loop for simulation and recording
        for k in range(time.size):
            while t < time[k]:

                # Update previous time and state
                t_old, x_old = t, x

                # Update current time and state
                u, x, is_jump = self.random_step(x, rng)
                t += u

                # Optional: record jump counts
                n_jumps[int(is_jump)] += 1

            # Record protein levels
            traj[k] = self.flow(time[k] - t_old, x_old)

        # Display info about jumps
        if verb:
            msg = (f'Exact simulation used {n_jumps.sum()} jumps '
                f'including {n_jumps[0]} phantom jumps '
                f'({100*n_jumps[0]/n_jumps.sum():.2f}%)')
            print(msg)

        return Trajectory(time, traj)


# Tests
if __name__ == '__main__':
    from models.networks import toggle_switch
    model = BurstyGRN(toggle_switch)
    # Simulation
    time = np.linspace(30, 50, 3)
    sim = model.simulate(time, verb=True, seed=0)
    print(sim.t)
    print(sim.x)
