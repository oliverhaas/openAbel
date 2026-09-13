# API reference

The public Python API of **openAbel** is a single class, `openAbel.Abel`. Everything else in the package is Cython
internals; the `.pxd` files are shipped, so the C-level functions can be `cimport`ed from other Cython modules, but
they are not a supported interface.

```python
import openAbel

openAbel.__version__  # e.g. "0.7.0"
openAbel.__all__  # ["Abel"]
```

## `openAbel.Abel`

```python
abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize, method=3, order=2, eps=1e1 * machineEpsilon)
```

Creates a transform plan for equispaced data of length `nData`. Creating the plan does the expensive preparation
(loading end-correction coefficients, building the FMM hierarchy); the plan is then reused for any number of
`execute` calls with the same grid.

| Parameter         | Type    | Description                                                                                                                                                                                                          |
| ----------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `nData`           | `int`   | Length of the data vector.                                                                                                                                                                                           |
| `forwardBackward` | `int`   | Which transform to perform: `-1` forward Abel transform, `1` backward (or inverse) Abel transform, `2` backward Abel transform with the derivative already supplied by the user, `-2` modified forward Abel transform. See [transform types](../transform-types.md). |
| `shift`           | `float` | Shift of the first sample away from 0 in positive direction, in units of `stepSize`. Usually `0.0` or `0.5`; the end-correction methods support only these two values.                                                 |
| `stepSize`        | `float` | Step size (or grid spacing) between two data points.                                                                                                                                                                 |
| `method`          | `int`   | Transform method: `0` desingularized trapezoidal rule, `1` Hansen-Law, `2` trapezoidal rule with end corrections, `3` Fast Multipole Method with end corrections (default). See [transform methods](../transform-methods.md). |
| `order`           | `int`   | Order of the end corrections for methods `2` and `3` (`0 < order < 20`, default `2`); ignored by methods `0` and `1`.                                                                                                |
| `eps`             | `float` | Target accuracy of the FMM far-field approximation (method `3` only); it sets the number of Chebyshev interpolation nodes. Must be at least the machine epsilon; defaults to ten times the machine epsilon.               |

Raises `ValueError` if a parameter has a non-viable value (for example `nData < 2`, `order <= 0`, or too few data points
for the requested order) and `NotImplementedError` if the chosen method does not support the given parameters (for
example an unknown `method`, a `shift` other than `0.0` or `0.5` with methods `2` and `3`, or the modified forward
transform with the Hansen-Law method).

### `Abel.execute`

```python
dataOut = abelObj.execute(dataIn, leftBoundary=0, rightBoundary=0)
```

Performs the transform and returns a new `numpy.ndarray` of `float64` with `nData` elements. `dataIn` is a
one-dimensional array of `float64` (any object supporting the buffer protocol with a double item type works); it is
copied, never modified.

| Parameter       | Type            | Description                                                                                                                                                                                                                            |
| --------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `dataIn`        | `numpy.ndarray` | Data vector.                                                                                                                                                                                                                           |
| `leftBoundary`  | `int`           | How the start of the data is handled: `0` data only given inside the integration interval, `1` data has odd symmetry around zero, `2` data has even symmetry around zero, `3` data is given outside the domain as well.                 |
| `rightBoundary` | `int`           | Almost the same as `leftBoundary`, only for the end of the data. `1` and `2` are not supported here.                                                                                                                                    |

With boundary value `3` the input is longer than `nData`: on that side it also carries the `(order - 1) // 2` samples
outside the domain that the end-correction stencil reaches into, plus `(orderFilter - 1) // 2` samples with
`orderFilter = order + 1 + order % 2` for the backward transform with numerical derivative (`forwardBackward=1`).
Samples beyond that are ignored; a shorter input raises `ValueError`. Method `0` behaves like `order = 1` here (one
extra sample per side for `forwardBackward=1`, none otherwise). For the Hansen-Law method (`method=1`) the boundary
arguments are ignored.

Raises `ValueError` for non-viable input and `NotImplementedError` for unsupported boundary combinations.

### Example

```python
import numpy as np

import openAbel

nData = 200
stepSize = 3.5 / (nData - 1)
x = np.arange(nData) * stepSize

abelObj = openAbel.Abel(nData, -1, 0.0, stepSize)  # forward transform, FMM with 2nd order end corrections
dataOut = abelObj.execute(np.exp(-(x**2)))
```

With the end-correction methods (`2` and `3`) the last sample of the result is exactly `0.0`: it is the truncated
transform at the truncation radius \(R\), where the integration interval has zero length. The
[examples](../examples/index.md) show how the truncation error behaves.
