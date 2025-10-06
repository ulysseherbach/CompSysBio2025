"""Various utility functions."""
import numpy as np


def check_time_points(time):
    """Check and return time points for trajectory simulations."""
    time = np.array(time, dtype=float, ndmin=1)
    if np.ndim(time) > 1:
        msg = 'Time points should either be scalar or 1D array.'
        raise ValueError(msg)
    if np.any(time != np.sort(time)):
        msg = 'Time points must be given in increasing order.'
        raise ValueError(msg)
    if np.any(time < 0):
        msg = 'Time points must be nonnegative.'
        raise ValueError(msg)
    return time


def check_init_state(init_state, shape=None):
    """Check and return initial state for trajectory simulations."""
    init_state = np.array(init_state, dtype=float)
    if np.any(init_state < 0):
        msg = 'Initial state must be nonnegative.'
        raise ValueError(msg)
    # Check for particular shape
    if (shape is not None) and (np.shape(init_state) != shape):
        if shape == ():
            msg = 'Initial state must be a scalar value.'
            raise ValueError(msg)
        msg = f'Initial state must have shape {shape}.'
        raise ValueError(msg)
    return init_state


# Tests
if __name__ == '__main__':
    time = np.linspace(0, 10, 11)
    check_time_points(time)
    state = 0, 1
    check_init_state(state, shape=(2,))
