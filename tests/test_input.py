import numpy as np
import pytest

import openabel


@pytest.mark.parametrize("method", [-1, 5])
def test_unknown_method_raises_not_implemented(method):
    with pytest.raises(NotImplementedError):
        openabel.Abel(10, -1, 0.0, 1.0, method=method)


@pytest.mark.parametrize("method", [2, 3])
def test_zero_order_raises_value_error(method):
    with pytest.raises(ValueError):
        openabel.Abel(10, -1, 0.0, 1.0, method=method, order=0)


@pytest.mark.parametrize("forward_backward", [-1, 1, 2, -2])
@pytest.mark.parametrize("method", [2, 3])
def test_unsupported_shift_raises_not_implemented(forward_backward, method):
    # Regression (method 3): the cleanup after this error freed uninitialised pointers and crashed the process.
    with pytest.raises(NotImplementedError):
        openabel.Abel(200, forward_backward, 0.25, 0.01, method=method)


@pytest.mark.parametrize("n_data", [-1, 0, 1])
@pytest.mark.parametrize("method", [0, 1, 2, 3])
def test_fewer_than_two_grid_points_raises_value_error(n_data, method):
    # Regression: n_data=0 aborted the process with method 1 and n_data=-1 reserved gigabytes.
    with pytest.raises(ValueError):
        openabel.Abel(n_data, -1, 0.0, 0.01, method=method)


@pytest.mark.parametrize(("n_data", "method", "order"), [(2, 0, 2), (5, 2, 4), (5, 3, 4)])
def test_too_few_data_points_for_method_raises_value_error(n_data, method, order):
    with pytest.raises(ValueError):
        openabel.Abel(n_data, -1, 0.0, 0.01, method=method, order=order)


def test_eps_below_machine_epsilon_raises_value_error():
    with pytest.raises(ValueError):
        openabel.Abel(200, -1, 0.0, 0.01, method=3, eps=0.0)


# Columns: forward_backward, method, order, n_outside (samples per side that boundary value 3 needs beyond n_data).
OUTSIDE_SAMPLES = (
    (-1, 0, 2, 0),
    (1, 0, 2, 1),
    (-1, 1, 2, 0),
    (1, 1, 2, 0),
    (-1, 2, 3, 1),
    (1, 2, 3, 3),
    (2, 2, 3, 1),
    (-1, 3, 3, 1),
    (1, 3, 3, 3),
    (1, 3, 2, 1),
)


@pytest.mark.parametrize(("forward_backward", "method", "order", "n_outside"), OUTSIDE_SAMPLES)
def test_outside_boundaries_short_input_raises_value_error(forward_backward, method, order, n_outside):
    # Regression: a short input crashed the process, and exactly n_data samples with boundary 3 read past the input.
    abel_obj = openabel.Abel(200, forward_backward, 0.0, 0.01, method=method, order=order)
    with pytest.raises(ValueError):
        abel_obj.execute(np.zeros(200 + 2 * n_outside - 1), left_boundary=3, right_boundary=3)


@pytest.mark.parametrize(("forward_backward", "method", "order", "n_outside"), OUTSIDE_SAMPLES)
def test_outside_boundaries_exact_input_returns_n_data_samples(forward_backward, method, order, n_outside):
    abel_obj = openabel.Abel(200, forward_backward, 0.0, 0.01, method=method, order=order)
    data_out = abel_obj.execute(np.zeros(200 + 2 * n_outside), left_boundary=3, right_boundary=3)
    assert data_out.shape == (200,)


@pytest.mark.parametrize(("forward_backward", "method", "order", "n_outside"), OUTSIDE_SAMPLES)
def test_left_outside_boundary_exact_input_returns_n_data_samples(forward_backward, method, order, n_outside):
    abel_obj = openabel.Abel(200, forward_backward, 0.0, 0.01, method=method, order=order)
    data_out = abel_obj.execute(np.zeros(200 + n_outside), left_boundary=3, right_boundary=0)
    assert data_out.shape == (200,)


@pytest.mark.parametrize("method", [0, 1, 2, 3])
def test_empty_input_raises_value_error(method):
    # Regression: an empty input corrupted the heap.
    abel_obj = openabel.Abel(200, -1, 0.0, 0.01, method=method)
    with pytest.raises(ValueError):
        abel_obj.execute(np.zeros(0))


def test_integer_input_raises_value_error():
    abel_obj = openabel.Abel(200, -1, 0.0, 0.01)
    with pytest.raises(ValueError):
        abel_obj.execute(np.zeros(200, dtype=np.int64))


def test_two_dimensional_input_raises_value_error():
    abel_obj = openabel.Abel(200, -1, 0.0, 0.01)
    with pytest.raises(ValueError):
        abel_obj.execute(np.zeros((200, 1)))


@pytest.mark.parametrize("method", [0, 2, 3])
def test_invalid_left_boundary_raises_not_implemented(method):
    abel_obj = openabel.Abel(200, -1, 0.0, 0.01, method=method)
    with pytest.raises(NotImplementedError):
        abel_obj.execute(np.zeros(200), left_boundary=4)


@pytest.mark.parametrize("right_boundary", [1, 2, 4])
@pytest.mark.parametrize("method", [0, 2, 3])
def test_unsupported_right_boundary_raises_not_implemented(right_boundary, method):
    abel_obj = openabel.Abel(200, -1, 0.0, 0.01, method=method)
    with pytest.raises(NotImplementedError):
        abel_obj.execute(np.zeros(200), right_boundary=right_boundary)


@pytest.mark.parametrize("method", [2, 3])
def test_order_without_coefficients_raises_value_error(method):
    with pytest.raises(ValueError):
        openabel.Abel(200, -1, 0.0, 0.01, method=method, order=20)
