"""Analytic Gaussian test pairs for every openAbel transform type on the truncated domain [0, R]."""

import numpy as np
from scipy import integrate, special

# Test function f(x) = exp(-x^2) on the grid x_i = (i + shift) * STEP_SIZE. Every reference is the transform of the
# truncated integral with R = x[-1], which is what openAbel computes; the last sample is 0 by construction.
N_DATA = 200
X_MAX = 3.5
STEP_SIZE = X_MAX / (N_DATA - 1)


def grid(*, shift: float) -> np.ndarray:
    return (np.arange(N_DATA) + shift) * STEP_SIZE


def relative_error(*, data_out: np.ndarray, reference: np.ndarray) -> float:
    """max |data_out - reference| / max |reference| over all but the last sample."""
    return float(np.max(np.abs(data_out[:-1] - reference[:-1])) / np.max(np.abs(reference)))


def _modified_forward_tail(*, y: float, R: float) -> float:
    """2 y^2 int_R^inf exp(-r^2) / (r^2 sqrt(r^2 - y^2)) dr, the part of the modified forward transform beyond R."""
    if y == 0.0:
        return 0.0
    value, _ = integrate.quad(lambda r: np.exp(-(r**2)) / (r**2 * np.sqrt(r**2 - y**2)), R, np.inf, limit=200)
    return 2.0 * y**2 * value


def input_samples(*, forward_backward: int, x: np.ndarray) -> np.ndarray:
    """The input function of the transform type ``forward_backward`` sampled at ``x`` (any points, also outside [0, R])."""
    g = np.exp(-(x**2))
    if forward_backward in (-1, -2):
        # forward and modified forward transform f(r) = exp(-r^2)
        return g
    if forward_backward == 1:
        # backward transform of F(y) = sqrt(pi) exp(-y^2); openAbel differentiates the input itself
        return np.sqrt(np.pi) * g
    if forward_backward == 2:
        # backward transform with the derivative F'(y) = -2 y sqrt(pi) exp(-y^2) supplied as input
        return -2.0 * x * np.sqrt(np.pi) * g
    msg = f"No analytic input for forward_backward={forward_backward}"
    raise ValueError(msg)


def analytic_pair(*, forward_backward: int, shift: float) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(data_in, expected data_out)`` for the transform type ``forward_backward`` on ``grid(shift)``."""
    x = grid(shift=shift)
    R = x[-1]
    g = np.exp(-(x**2))
    data_in = input_samples(forward_backward=forward_backward, x=x)
    truncated_erf = special.erf(np.sqrt(np.maximum(R**2 - x**2, 0.0)))
    if forward_backward == -1:
        # forward: F(y) = 2 int_y^R r exp(-r^2) / sqrt(r^2 - y^2) dr = sqrt(pi) exp(-y^2) erf(sqrt(R^2 - y^2))
        return data_in, np.sqrt(np.pi) * g * truncated_erf
    if forward_backward in (1, 2):
        # backward: f(r) = exp(-r^2), truncated at R
        return data_in, g * truncated_erf
    # modified forward: H(y) = 2 int_y^inf exp(-r^2) y^2 / (r^2 sqrt(r^2 - y^2)) dr
    #                        = y^2 exp(-y^2) [K1e(y^2/2) - K0e(y^2/2)], H(0) = 2, minus the tail beyond R
    h = np.full_like(x, 2.0)
    xp = x[x > 0.0]
    h[x > 0.0] = xp**2 * np.exp(-(xp**2)) * (special.k1e(0.5 * xp**2) - special.k0e(0.5 * xp**2))
    tail = np.array([*(_modified_forward_tail(y=y, R=R) for y in x[:-1]), h[-1]])
    return data_in, h - tail
