import ctypes
from types import MappingProxyType

import numpy as np
import pytest
from analytic import N_DATA, STEP_SIZE, analytic_pair, input_samples, relative_error

import openabel

# Relative-error tolerances: the larger of the two shifts' errors measured on the Cython 3 build of 2026-09-13, times
# 5, rounded up to the next power of ten (spec section 3). Exceptions: order 10 sits 10-100x above the rule because
# the measured errors (1e-15 to 7e-13) are within reach of BLAS noise; (fb=2, method=1) and (fb=1, method=0) are 3x
# the measured 6.4e-2 and 6.8e-2 because the rule would give 1.
# Key: (forward_backward, order) for the end-correction methods 2 and 3.
END_CORRECTION_TOLERANCE = MappingProxyType(
    {
        (-1, 1): 1e-2,
        (-1, 2): 1e-4,
        (-1, 3): 1e-6,
        (-1, 5): 1e-9,
        (-1, 10): 1e-13,
        (1, 1): 1e-1,
        (1, 2): 1e-1,  # shift 0.5 measures 1.14e-2, shift 0 measures 2.2e-4
        (1, 3): 1e-4,
        (1, 5): 1e-6,
        (1, 10): 1e-11,
        (2, 1): 1e-1,
        (2, 2): 1e-4,
        (2, 3): 1e-4,
        (2, 5): 1e-7,
        (2, 10): 1e-13,
        (-2, 1): 1e-2,
        (-2, 2): 1e-3,
        (-2, 3): 1e-6,
        (-2, 5): 1e-9,
        (-2, 10): 1e-11,
    },
)
# Key: (forward_backward, method) for the single-order methods 0 (desingularised trapezoidal) and 1 (Hansen-Law).
SINGLE_ORDER_TOLERANCE = MappingProxyType(
    {
        (-1, 0): 1e-2,
        (-1, 1): 1e-2,
        # measured 6.8e-2 at shift 0, 9.5e-3 at shift 0.5; regression: raised FileNotFoundError before 0.7.0
        (1, 0): 2e-1,
        (1, 1): 1e-1,
        (2, 0): 1e-1,
        (2, 1): 2e-1,  # measured 6.4e-2
        (-2, 0): 1e-2,
    },
)


@pytest.mark.parametrize("method", [2, 3])
@pytest.mark.parametrize("shift", [0.0, 0.5])
@pytest.mark.parametrize(
    ("forward_backward", "order", "tolerance"),
    [(forward_backward, order, tolerance) for (forward_backward, order), tolerance in END_CORRECTION_TOLERANCE.items()],
)
def test_end_correction_methods_match_analytic_transform(forward_backward, order, tolerance, shift, method):
    data_in, reference = analytic_pair(forward_backward=forward_backward, shift=shift)
    data_out = openabel.Abel(N_DATA, forward_backward, shift, STEP_SIZE, method=method, order=order).execute(data_in)
    assert data_out[-1] == 0.0
    assert relative_error(data_out=data_out, reference=reference) < tolerance


@pytest.mark.parametrize("shift", [0.0, 0.5])
@pytest.mark.parametrize(
    ("forward_backward", "method", "tolerance"),
    [(forward_backward, method, tolerance) for (forward_backward, method), tolerance in SINGLE_ORDER_TOLERANCE.items()],
)
def test_single_order_methods_match_analytic_transform(forward_backward, method, tolerance, shift):
    data_in, reference = analytic_pair(forward_backward=forward_backward, shift=shift)
    data_out = openabel.Abel(N_DATA, forward_backward, shift, STEP_SIZE, method=method).execute(data_in)
    assert data_out[-1] == 0.0
    assert np.isfinite(data_out).all()
    assert relative_error(data_out=data_out, reference=reference) < tolerance


@pytest.mark.parametrize("shift", [0.25, 1.0, 3.0])
@pytest.mark.parametrize(("forward_backward", "tolerance"), [(-1, 1e-2), (1, 5e-2), (2, 2e-2), (-2, 1e-2)])
def test_desingularised_trapezoidal_supports_any_positive_shift(forward_backward, tolerance, shift):
    data_in, reference = analytic_pair(forward_backward=forward_backward, shift=shift)
    data_out = openabel.Abel(N_DATA, forward_backward, shift, STEP_SIZE, method=0).execute(data_in)
    assert np.isfinite(data_out).all()
    assert relative_error(data_out=data_out, reference=reference) < tolerance


def outside_samples_per_side(*, forward_backward: int, order: int) -> int:
    """Samples outside the domain that boundary value 3 consumes per side: the half widths of the end-correction
    stencil and, for the backward transform with numerical derivative, of the derivative filter."""
    order_filter = order + 1 + order % 2 if forward_backward == 1 else 1
    return (order - 1) // 2 + (order_filter - 1) // 2


def test_fmm_backward_order2_ignores_input_beyond_stencil_reach():
    # Regression: even orders read one sample past the stencil reach and the FMM multiplied it (or, with boundary 0,
    # an uninitialised buffer element) by a zero coefficient, so a NaN there poisoned the result.
    n_outside = outside_samples_per_side(forward_backward=1, order=2)
    x = np.arange(-n_outside, N_DATA + n_outside + 1) * STEP_SIZE
    data_in = input_samples(forward_backward=1, x=x)
    data_in[-1] = np.nan  # one sample beyond what boundary value 3 needs
    abel_obj = openabel.Abel(N_DATA, 1, 0.0, STEP_SIZE, method=3, order=2)
    data_out = abel_obj.execute(data_in, left_boundary=3, right_boundary=3)
    assert np.isfinite(data_out).all()


@pytest.mark.parametrize("method", [2, 3])
@pytest.mark.parametrize("shift", [0.0, 0.5])
@pytest.mark.parametrize(
    ("forward_backward", "order", "tolerance"),
    [(forward_backward, order, tolerance) for (forward_backward, order), tolerance in END_CORRECTION_TOLERANCE.items()],
)
def test_end_correction_methods_outside_samples_match_analytic_transform(
    forward_backward,
    order,
    tolerance,
    shift,
    method,
):
    n_outside = outside_samples_per_side(forward_backward=forward_backward, order=order)
    x = (np.arange(-n_outside, N_DATA + n_outside) + shift) * STEP_SIZE
    data_in = input_samples(forward_backward=forward_backward, x=x)
    _, reference = analytic_pair(forward_backward=forward_backward, shift=shift)
    abel_obj = openabel.Abel(N_DATA, forward_backward, shift, STEP_SIZE, method=method, order=order)
    data_out = abel_obj.execute(data_in, left_boundary=3, right_boundary=3)
    assert data_out.shape == (N_DATA,)
    assert data_out[-1] == 0.0
    assert relative_error(data_out=data_out, reference=reference) < tolerance


def symmetric_left_boundary(*, forward_backward: int) -> int:
    """1 (odd) for the derivative input of forward_backward=2, 2 (even) for the other Gaussian inputs."""
    return 1 if forward_backward == 2 else 2


@pytest.mark.parametrize("method", [2, 3])
@pytest.mark.parametrize("shift", [0.0, 0.5])
@pytest.mark.parametrize(
    ("forward_backward", "order", "tolerance"),
    [(forward_backward, order, tolerance) for (forward_backward, order), tolerance in END_CORRECTION_TOLERANCE.items()],
)
def test_end_correction_methods_symmetric_boundary_match_analytic_transform(
    forward_backward,
    order,
    tolerance,
    shift,
    method,
):
    data_in, reference = analytic_pair(forward_backward=forward_backward, shift=shift)
    abel_obj = openabel.Abel(N_DATA, forward_backward, shift, STEP_SIZE, method=method, order=order)
    data_out = abel_obj.execute(data_in, left_boundary=symmetric_left_boundary(forward_backward=forward_backward))
    assert data_out[-1] == 0.0
    assert relative_error(data_out=data_out, reference=reference) < tolerance


@pytest.mark.parametrize("shift", [0.0, 0.5])
def test_desingularised_trapezoidal_even_boundary_improves_numerical_derivative(shift):
    data_in, reference = analytic_pair(forward_backward=1, shift=shift)
    abel_obj = openabel.Abel(N_DATA, 1, shift, STEP_SIZE, method=0)
    extrapolated = relative_error(data_out=abel_obj.execute(data_in, left_boundary=0), reference=reference)
    mirrored = relative_error(data_out=abel_obj.execute(data_in, left_boundary=2), reference=reference)
    assert mirrored < extrapolated
    assert mirrored < SINGLE_ORDER_TOLERANCE[1, 0]


@pytest.mark.parametrize("eps", [1e-3, 1e-6, 1e-10])
@pytest.mark.parametrize("forward_backward", [-1, 1, 2, -2])
def test_fmm_eps_bounds_the_deviation_from_the_default_eps(forward_backward, eps):
    data_in, _ = analytic_pair(forward_backward=forward_backward, shift=0.0)
    default = openabel.Abel(N_DATA, forward_backward, 0.0, STEP_SIZE).execute(data_in)
    loose = openabel.Abel(N_DATA, forward_backward, 0.0, STEP_SIZE, eps=eps).execute(data_in)
    deviation = np.max(np.abs(loose - default)) / np.max(np.abs(default))
    assert 0.0 < deviation < eps


@pytest.mark.parametrize("n_data", [4, 5, 6, 7, 8])
def test_fmm_with_few_data_points_matches_trapezoidal_without_blas_complaints(n_data, capfd):
    x = np.arange(n_data) * STEP_SIZE
    data_in = input_samples(forward_backward=-1, x=x)
    trapezoidal = openabel.Abel(n_data, -1, 0.0, STEP_SIZE, method=2).execute(data_in)
    fmm = openabel.Abel(n_data, -1, 0.0, STEP_SIZE, method=3).execute(data_in)
    ctypes.CDLL(None).fflush(None)
    captured = capfd.readouterr()
    assert captured.out == ""
    assert captured.err == ""
    np.testing.assert_allclose(fmm, trapezoidal, rtol=1e-8, atol=1e-8)


@pytest.mark.parametrize("order", [1, 2, 3, 5, 10])
def test_modified_forward_half_shift_trapezoidal_and_fmm_agree(order):
    # Regression: the half-shift coefficient key had a trailing underscore; method 2 raised KeyError, method 3 crashed.
    data_in, _ = analytic_pair(forward_backward=-2, shift=0.5)
    trapezoidal = openabel.Abel(N_DATA, -2, 0.5, STEP_SIZE, method=2, order=order).execute(data_in)
    fmm = openabel.Abel(N_DATA, -2, 0.5, STEP_SIZE, method=3, order=order).execute(data_in)
    np.testing.assert_allclose(fmm, trapezoidal, rtol=1e-8, atol=1e-8)


@pytest.mark.parametrize("shift", [0.0, 0.5])
def test_hansen_law_modified_forward_raises_not_implemented(shift):
    # Regression: Cython 0.29 swallowed this NotImplementedError (no except clause) and returned a copy of the input.
    with pytest.raises(NotImplementedError):
        openabel.Abel(N_DATA, -2, shift, STEP_SIZE, method=1)


def test_hansen_law_works_after_failed_construction():
    # Regression: the failed construction above leaked its plan data; a following transform must be unaffected.
    with pytest.raises(NotImplementedError):
        openabel.Abel(N_DATA, -2, 0.0, STEP_SIZE, method=1)
    data_in, reference = analytic_pair(forward_backward=-1, shift=0.0)
    data_out = openabel.Abel(N_DATA, -1, 0.0, STEP_SIZE, method=1).execute(data_in)
    assert relative_error(data_out=data_out, reference=reference) < SINGLE_ORDER_TOLERANCE[-1, 1]
