from types import MappingProxyType

import numpy as np
import pytest
from analytic import N_DATA, STEP_SIZE, analyticPair, relativeError

import openAbel

# Relative-error tolerances: the larger of the two shifts' errors measured on the Cython 3 build of 2026-09-13, times
# 5, rounded up to the next power of ten (spec section 3). Exceptions: order 10 sits 10-100x above the rule because
# the measured errors (1e-15 to 7e-13) are within reach of BLAS noise; (fb=2, method=1) and (fb=1, method=0) are 3x
# the measured 6.4e-2 and 6.8e-2 because the rule would give 1.
# Key: (forwardBackward, order) for the end-correction methods 2 and 3.
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
# Key: (forwardBackward, method) for the single-order methods 0 (desingularised trapezoidal) and 1 (Hansen-Law).
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
    ("forwardBackward", "order", "tolerance"),
    [(forwardBackward, order, tolerance) for (forwardBackward, order), tolerance in END_CORRECTION_TOLERANCE.items()],
)
def test_endCorrectionMethods_matchAnalyticTransform(forwardBackward, order, tolerance, shift, method):
    dataIn, reference = analyticPair(forwardBackward=forwardBackward, shift=shift)
    dataOut = openAbel.Abel(N_DATA, forwardBackward, shift, STEP_SIZE, method=method, order=order).execute(dataIn)
    assert dataOut[-1] == 0.0
    assert relativeError(dataOut=dataOut, reference=reference) < tolerance


@pytest.mark.parametrize("shift", [0.0, 0.5])
@pytest.mark.parametrize(
    ("forwardBackward", "method", "tolerance"),
    [(forwardBackward, method, tolerance) for (forwardBackward, method), tolerance in SINGLE_ORDER_TOLERANCE.items()],
)
def test_singleOrderMethods_matchAnalyticTransform(forwardBackward, method, tolerance, shift):
    dataIn, reference = analyticPair(forwardBackward=forwardBackward, shift=shift)
    dataOut = openAbel.Abel(N_DATA, forwardBackward, shift, STEP_SIZE, method=method).execute(dataIn)
    assert dataOut[-1] == 0.0
    assert np.isfinite(dataOut).all()
    assert relativeError(dataOut=dataOut, reference=reference) < tolerance


@pytest.mark.parametrize("order", [1, 2, 3, 5, 10])
def test_modifiedForwardHalfShift_trapezoidalAndFmmAgree(order):
    # Regression: the half-shift coefficient key had a trailing underscore; method 2 raised KeyError, method 3 crashed.
    dataIn, _ = analyticPair(forwardBackward=-2, shift=0.5)
    trapezoidal = openAbel.Abel(N_DATA, -2, 0.5, STEP_SIZE, method=2, order=order).execute(dataIn)
    fmm = openAbel.Abel(N_DATA, -2, 0.5, STEP_SIZE, method=3, order=order).execute(dataIn)
    np.testing.assert_allclose(fmm, trapezoidal, rtol=1e-8, atol=1e-8)


@pytest.mark.parametrize("shift", [0.0, 0.5])
def test_hansenLawModifiedForward_raisesNotImplemented(shift):
    # Regression: Cython 0.29 swallowed this NotImplementedError (no except clause) and returned a copy of the input.
    with pytest.raises(NotImplementedError):
        openAbel.Abel(N_DATA, -2, shift, STEP_SIZE, method=1)


def test_hansenLaw_worksAfterFailedConstruction():
    # Regression: the failed construction above leaked its plan data; a following transform must be unaffected.
    with pytest.raises(NotImplementedError):
        openAbel.Abel(N_DATA, -2, 0.0, STEP_SIZE, method=1)
    dataIn, reference = analyticPair(forwardBackward=-1, shift=0.0)
    dataOut = openAbel.Abel(N_DATA, -1, 0.0, STEP_SIZE, method=1).execute(dataIn)
    assert relativeError(dataOut=dataOut, reference=reference) < SINGLE_ORDER_TOLERANCE[-1, 1]
