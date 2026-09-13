# API reference

The public Python API of **openAbel** is a single class, `openabel.Abel`. Everything else in the package is Cython
internals; the `.pxd` files are shipped, so the C-level functions can be `cimport`ed from other Cython modules, but
they are not a supported interface.

```python
import openabel

openabel.__version__  # e.g. "0.7.0"
openabel.__all__  # ["Abel"]
```

## `openabel.Abel`

```python
abel_obj = openabel.Abel(n_data, forward_backward, shift, step_size, method=3, order=2, eps=1e1 * machine_epsilon)
```

Creates a transform plan for equispaced data of length `n_data`. Creating the plan does the expensive preparation
(loading end-correction coefficients, building the FMM hierarchy); the plan is then reused for any number of
`execute` calls with the same grid.

| Parameter         | Type    | Description                                                                                                                                                                                                          |
| ----------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `n_data`           | `int`   | Length of the data vector.                                                                                                                                                                                           |
| `forward_backward` | `int`   | Which transform to perform: `-1` forward Abel transform, `1` backward (or inverse) Abel transform, `2` backward Abel transform with the derivative already supplied by the user, `-2` modified forward Abel transform. See [transform types](../transform-types.md). |
| `shift`           | `float` | Shift of the first sample away from 0 in positive direction, in units of `step_size`. Usually `0.0` or `0.5`; the end-correction methods support only these two values.                                                 |
| `step_size`        | `float` | Step size (or grid spacing) between two data points.                                                                                                                                                                 |
| `method`          | `int`   | Transform method: `0` desingularized trapezoidal rule, `1` Hansen-Law, `2` trapezoidal rule with end corrections, `3` Fast Multipole Method with end corrections (default). See [transform methods](../transform-methods.md). |
| `order`           | `int`   | Order of the end corrections for methods `2` and `3` (`0 < order < 20`, default `2`); ignored by methods `0` and `1`.                                                                                                |
| `eps`             | `float` | Target accuracy of the FMM far-field approximation (method `3` only); it sets the number of Chebyshev interpolation nodes. Must be at least the machine epsilon; defaults to ten times the machine epsilon.               |

Raises `ValueError` if a parameter has a non-viable value (for example `n_data < 2`, `order <= 0`, an `order` without
coefficient tables (`order >= 20`), or too few data points for the requested order) and `NotImplementedError` if the
chosen method does not support the given parameters (for example an unknown `method`, a `shift` other than `0.0` or
`0.5` with methods `2` and `3`, or the modified forward transform with the Hansen-Law method).

### `Abel.execute`

```python
data_out = abel_obj.execute(data_in, left_boundary=0, right_boundary=0)
```

Performs the transform and returns a new `numpy.ndarray` of `float64` with `n_data` elements. `data_in` is a
one-dimensional array of `float64` (any object supporting the buffer protocol with a double item type works); it is
copied, never modified.

| Parameter       | Type            | Description                                                                                                                                                                                                                            |
| --------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `data_in`        | `numpy.ndarray` | Data vector.                                                                                                                                                                                                                           |
| `left_boundary`  | `int`           | How the start of the data is handled: `0` data only given inside the integration interval, `1` data has odd symmetry around zero, `2` data has even symmetry around zero, `3` data is given outside the domain as well.                 |
| `right_boundary` | `int`           | Almost the same as `left_boundary`, only for the end of the data. `1` and `2` are not supported here.                                                                                                                                    |

With boundary value `3` the input is longer than `n_data`: on that side it also carries the `(order - 1) // 2` samples
outside the domain that the end-correction stencil reaches into, plus `(order_filter - 1) // 2` samples with
`order_filter = order + 1 + order % 2` for the backward transform with numerical derivative (`forward_backward=1`).
Samples beyond that are ignored; a shorter input raises `ValueError`. Method `0` behaves like `order = 1` here (one
extra sample per side for `forward_backward=1`, none otherwise). For the Hansen-Law method (`method=1`) the boundary
arguments are ignored.

Raises `ValueError` for non-viable input and `NotImplementedError` for unsupported boundary combinations.

### Example

```python
import numpy as np

import openabel

n_data = 200
step_size = 3.5 / (n_data - 1)
x = np.arange(n_data) * step_size

abel_obj = openabel.Abel(n_data, -1, 0.0, step_size)  # forward transform, FMM with 2nd order end corrections
data_out = abel_obj.execute(np.exp(-(x**2)))
```

With the end-correction methods (`2` and `3`) the last sample of the result is exactly `0.0`: it is the truncated
transform at the truncation radius \(R\), where the integration interval has zero length. The
[examples](../examples/index.md) show how the truncation error behaves.
