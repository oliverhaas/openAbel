import pytest

import openAbel


@pytest.mark.parametrize("method", [-1, 5])
def test_unknownMethod_raisesNotImplemented(method):
    with pytest.raises(NotImplementedError):
        openAbel.Abel(10, -1, 0.0, 1.0, method=method)


@pytest.mark.parametrize("method", [2, 3])
def test_zeroOrder_raisesValueError(method):
    with pytest.raises(ValueError):
        openAbel.Abel(10, -1, 0.0, 1.0, method=method, order=0)


@pytest.mark.parametrize("forwardBackward", [-1, 1, 2, -2])
@pytest.mark.parametrize("method", [2, 3])
def test_unsupportedShift_raisesNotImplemented(forwardBackward, method):
    # Regression (method 3): md.ltp and md.direct0 were not NULL-initialised, so the cleanup after this error freed
    # garbage and the process crashed instead of raising.
    with pytest.raises(NotImplementedError):
        openAbel.Abel(200, forwardBackward, 0.25, 0.01, method=method)
