"""Simulation of the Bursty model for a single gene with no feedback."""
import numpy as np
import models._utils as utils


class Trajectory:
    """Store single-cell gene expression trajectories."""

    def __init__(self, t, x):
        self.t = t  # Time points
        self.x = x  # Protein levels


class BurstyBase:
    """Bursty model for a single gene with no feedback."""

    def __init__(self,
        burst_size=1.0,
        burst_frequency=1.0,
        degradation_rate=1.0):

        # Set model parameters
        self.burst_size = burst_size
        self.burst_frequency = burst_frequency
        self.degradation_rate = degradation_rate

    def simulate(self, time, init_state=0.0, seed=None, verb=False):
        """Perform exact simulation (extracted at given time points)."""
        burst_size = self.burst_size
        burst_frequency = self.burst_frequency
        degradation_rate = self.degradation_rate

        # Check simulation parameters
        time = utils.check_time_points(time).reshape((-1,))
        init_state = utils.check_init_state(init_state, shape=())

        # Define a random generator
        rng = np.random.default_rng(seed)

        # Sample total number of bursts
        n = rng.poisson(burst_frequency * time[-1])

        # Sample burst times and heights
        t = rng.uniform(low=0, high=time[-1], size=n)
        h = rng.exponential(scale=burst_size, size=n)

        # Build the post-jump embedded Markov chain
        t = np.append(0, np.sort(t))
        x = np.zeros(n + 1)
        x[0] = init_state
        for k in range(n):
            dt = t[k+1] - t[k]
            x[k+1] = x[k] * np.exp(- degradation_rate * dt) + h[k]

        # Extract states at user time points
        traj = np.zeros(time.size)
        for i, u in enumerate(time):
            k = np.argwhere(t <= u)[-1][0]
            dt = u - t[k]
            traj[i] = x[k] * np.exp(- degradation_rate * dt)
        if verb:
            print(f'Simulation generated {n} jumps.')

        return Trajectory(time, traj)


# Tests
if __name__ == '__main__':
    model = BurstyBase()
    # Simulation
    time = np.linspace(0, 50, 6)
    sim = model.simulate(time, verb=True)
    print(sim.t)
    print(sim.x)
