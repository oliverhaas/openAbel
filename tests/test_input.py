import numpy as np
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


@pytest.mark.parametrize("nData", [-1, 0, 1])
@pytest.mark.parametrize("method", [0, 1, 2, 3])
def test_fewerThanTwoGridPoints_raisesValueError(nData, method):
    # Regression: nData=0 with method 1 wrote past a zero-size block and aborted the process, nData=-1 reserved
    # gigabytes.
    with pytest.raises(ValueError):
        openAbel.Abel(nData, -1, 0.0, 0.01, method=method)


@pytest.mark.parametrize(("nData", "method", "order"), [(2, 0, 2), (5, 2, 4), (5, 3, 4)])
def test_tooFewDataPointsForMethod_raisesValueError(nData, method, order):
    with pytest.raises(ValueError):
        openAbel.Abel(nData, -1, 0.0, 0.01, method=method, order=order)


def test_epsBelowMachineEpsilon_raisesValueError():
    with pytest.raises(ValueError):
        openAbel.Abel(200, -1, 0.0, 0.01, method=3, eps=0.0)


# Samples per side that boundary value 3 needs beyond nData: the half widths of the end-correction stencil and, for
# the backward transform with numerical derivative, of the derivative filter. Method 0 has a first-order stencil,
# method 1 ignores the boundary values. Columns: forwardBackward, method, order, nOutside.
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


@pytest.mark.parametrize(("forwardBackward", "method", "order", "nOutside"), OUTSIDE_SAMPLES)
def test_outsideBoundaries_shortInput_raisesValueError(forwardBackward, method, order, nOutside):
    # Regression: a short input crashed the process, and boundary value 3 with exactly nData samples read past the
    # input and returned garbage.
    abelObj = openAbel.Abel(200, forwardBackward, 0.0, 0.01, method=method, order=order)
    with pytest.raises(ValueError):
        abelObj.execute(np.zeros(200 + 2 * nOutside - 1), leftBoundary=3, rightBoundary=3)


@pytest.mark.parametrize(("forwardBackward", "method", "order", "nOutside"), OUTSIDE_SAMPLES)
def test_outsideBoundaries_exactInput_returnsNDataSamples(forwardBackward, method, order, nOutside):
    abelObj = openAbel.Abel(200, forwardBackward, 0.0, 0.01, method=method, order=order)
    dataOut = abelObj.execute(np.zeros(200 + 2 * nOutside), leftBoundary=3, rightBoundary=3)
    assert dataOut.shape == (200,)


@pytest.mark.parametrize(("forwardBackward", "method", "order", "nOutside"), OUTSIDE_SAMPLES)
def test_leftOutsideBoundary_exactInput_returnsNDataSamples(forwardBackward, method, order, nOutside):
    abelObj = openAbel.Abel(200, forwardBackward, 0.0, 0.01, method=method, order=order)
    dataOut = abelObj.execute(np.zeros(200 + nOutside), leftBoundary=3, rightBoundary=0)
    assert dataOut.shape == (200,)


@pytest.mark.parametrize("method", [0, 1, 2, 3])
def test_emptyInput_raisesValueError(method):
    # Regression: an empty input corrupted the heap.
    abelObj = openAbel.Abel(200, -1, 0.0, 0.01, method=method)
    with pytest.raises(ValueError):
        abelObj.execute(np.zeros(0))


def test_integerInput_raisesValueError():
    abelObj = openAbel.Abel(200, -1, 0.0, 0.01)
    with pytest.raises(ValueError):
        abelObj.execute(np.zeros(200, dtype=np.int64))


def test_twoDimensionalInput_raisesValueError():
    abelObj = openAbel.Abel(200, -1, 0.0, 0.01)
    with pytest.raises(ValueError):
        abelObj.execute(np.zeros((200, 1)))


@pytest.mark.parametrize("method", [0, 2, 3])
def test_invalidLeftBoundary_raisesNotImplemented(method):
    abelObj = openAbel.Abel(200, -1, 0.0, 0.01, method=method)
    with pytest.raises(NotImplementedError):
        abelObj.execute(np.zeros(200), leftBoundary=4)


@pytest.mark.parametrize("rightBoundary", [1, 2, 4])
@pytest.mark.parametrize("method", [0, 2, 3])
def test_unsupportedRightBoundary_raisesNotImplemented(rightBoundary, method):
    abelObj = openAbel.Abel(200, -1, 0.0, 0.01, method=method)
    with pytest.raises(NotImplementedError):
        abelObj.execute(np.zeros(200), rightBoundary=rightBoundary)
