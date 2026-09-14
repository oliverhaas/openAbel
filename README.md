# openAbel

[![CI](https://github.com/oliverhaas/openAbel/actions/workflows/ci.yml/badge.svg)](https://github.com/oliverhaas/openAbel/actions/workflows/ci.yml)

Fast Abel transforms of equispaced data in Python, with all calculations done in Cython.

## Introduction

The main goal of **openAbel** is to provide fast and efficient Abel transforms of equispaced data
in Python with all actual calculations done in Cython. The most useful methods implemented in this
module for that purpose use the Fast Multipole Method combined with arbitrary order end correction
of the trapezoidal rule to achieve small errors and fast convergence, as well as linear computational
complexity. A couple of other methods are implemented for comparisons. The Abel transform can be
used from Python with numpy arrays or from Cython using pointers.

## Quick start

Requirements: Python >= 3.12 on Linux or macOS. numpy and scipy are installed automatically.

Once released on PyPI:

```bash
pip install openabel
```

Until then, install from the repository (this compiles the Cython extensions, so a C compiler is needed):

```bash
pip install git+https://github.com/oliverhaas/openAbel
```

A forward transform of a Gaussian sampled on 200 points, the first sample at `x = 0`:

```python
import numpy as np
import openabel

n_data = 200
step_size = 3.5 / (n_data - 1)
x = np.arange(n_data) * step_size

abel_obj = openabel.Abel(n_data, -1, 0.0, step_size)
data_out = abel_obj.execute(np.exp(-(x**2)))
```

The [examples](https://github.com/oliverhaas/openAbel/tree/main/examples) show the transform types and
methods in more detail, starting with `example000_simple_forward.py`. They need matplotlib (and
`example005` PyAbel), which a checkout provides as the `examples` dependency group:

```bash
uv run --group examples python examples/example000_simple_forward.py
```

## Development

```bash
uv sync                    # builds the extensions into .venv (editable install)
uv run pytest
uv run pre-commit install  # ruff and the other hooks on commit, ty on push
```

After editing a `.pyx` or `.pxd` file, rebuild with `uv sync --reinstall-package openabel`.

## Issues

If there are any issues, bugs or feature requests just let me know. As of now there are some gaps in
the implementation, e.g. not all transform types are available in all methods, but since the default
method vastly outperforms every other method anyway it's not really a pressing issue.

## Transform methods

For the default and most important method of **openAbel** we adapted the Chebyshev interpolation Fast
Multipole Method (FMM) as described by
[Tausch](https://link.springer.com/chapter/10.1007/978-3-642-25670-7_6) and calculated end corrections
specifically for the Abel transform similar to
[Kapur](https://epubs.siam.org/doi/abs/10.1137/S0036142995287847). If data points outside of the
integration interval can be provided these end corrections are arbitrary order stable and we provide
coefficients up to 19th order, otherwise it's recommended to use at most 5th order. The FMM leads to a
linear *O(N)* computational complexity algorithm.

In both error and computational complexity there is no better existing method for the intended purpose
to my knowledge. I should really stress that there are dozens of publications and methods out there
which claim to be fast and/or accurate, but don't get anywhere close to **openAbel** in those aspects.

For more information see the documentation on the
[transform methods](https://github.com/oliverhaas/openAbel/blob/main/docs/transform-methods.md) and the
[examples](https://github.com/oliverhaas/openAbel/blob/main/docs/examples/index.md).

## Copyright and license

Copyright 2016-2026 Oliver Sebastian Haas.

The code **openAbel** is published under the GNU GPL version 3. This program is free software; you can
redistribute it and/or modify it under the terms of the GNU General Public License as published by the
Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

For more information see the GNU General Public License copy provided in this repository:
[LICENSE](https://github.com/oliverhaas/openAbel/blob/main/LICENSE).
