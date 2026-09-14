# Changelog

## 0.7.0 (2026-09-14)

### Packaging and tooling

- Python >= 3.12; wheels for CPython 3.12, 3.13, 3.14 and the free-threaded 3.14t on Linux x86_64 and macOS arm64.
- `pyproject.toml` (PEP 621) with a `src/` layout; `setup.py` only compiles the Cython extensions. Version `0.7.0` is
  exposed as `openabel.__version__`.
- Cython 3, numpy 2 and scipy >= 1.13.
- pytest replaces nose; ruff, ty and pre-commit; GitHub Actions CI (Linux 3.12-3.14t, macOS 3.14) replaces Travis.
- Documentation ported from Sphinx to mkdocs-material.
- The sdist ships the tests and no generated C files; the `examples` dependency group installs matplotlib and PyAbel
  for the example scripts.

### Fixed

- Every `Abel(...)` construction died with exit status 255 on glibc >= 2.38 (e.g. Ubuntu 24.04): the internal
  allocator requested alignment 0 from `aligned_alloc` and then called `exit` on the NULL it got back. The allocator
  now uses 64-byte alignment, rounds the size up as C11 requires, and raises `MemoryError` on failure.
- Backward transform with `method=0` raised `FileNotFoundError` (wrong coefficient path).
- Modified forward transform (`forward_backward=-2`) with `shift=0.5` raised `KeyError` with `method=2` and crashed
  with `method=3` (wrong coefficient key).
- An unsupported `shift` with `method=3` crashed the process instead of raising `NotImplementedError` (uninitialised
  pointers were freed during cleanup).
- The Hansen-Law method (`method=1`) silently returned a copy of the input for the modified forward transform and
  leaked memory; it now raises `NotImplementedError`.
- With even `order` (including the default `2`) the end-correction methods (`2` and `3`) read one sample past the end
  of the input, and the FMM (`method=3`) multiplied an uninitialised buffer element by a zero coefficient. When that
  element happened to hold a NaN or Inf bit pattern the result was garbage, so the backward transform with `method=3`
  failed sporadically. The buffers are now sized exactly; results are unchanged otherwise.
- `method=0` leaked a small allocation per `Abel(...)` construction.
- The backward transforms with `method=0` and a `shift` other than `0.0` or `0.5` left the first desingularisation
  weight uninitialised, so the first output sample was off by up to 30% (or garbage).
- `Abel(...)` accepted `n_data < 2` and `execute` accepted inputs shorter than the plan needs; depending on the
  method the process crashed or the result was garbage. Both now raise `ValueError`, and the message names the
  required length (the boundary value `3` rule in the API reference).
- An invalid `left_boundary` value leaked two temporary buffers with methods `0` and `2`.
- An `order` without coefficient tables (20 and above) raised `KeyError`; it now raises `ValueError`.
- A `step_size` that is not positive or a negative `shift` silently produced NaN, Inf or mirrored results; both now
  raise `ValueError`.
- The FMM (`method=3`) with fewer than eight data points called DGEMM with a negative dimension; the results were
  right, but OpenBLAS printed an "illegal value" complaint.
- The FMM (`method=3`) sized its coefficient tables and indexed them with `int` arithmetic, which overflows for
  `n_data` above roughly 1e7; the sizes and offsets are now computed in `size_t`.

### Changed

- The package is imported and installed as `openabel` (`pip install openabel`, `import openabel`; PyPI treats it as
  the same project as `openAbel`), and every name follows PEP 8 snake_case: `Abel(n_data, forward_backward, shift,
  step_size, method=3, order=2, eps=...)` and `execute(data_in, left_boundary=0, right_boundary=0)`. Positional
  calls from 0.6 keep working once the import is updated; keyword arguments need the new names. The Cython module
  `hansen_law`, the internals and the coefficient data files are renamed the same way; the unused `math_fun` module
  and the unused `coeffs_inv_sqrt_sing` coefficient tables are gone.
- Cython 3 build: exception clauses moved after `nogil`, `cpow=True` keeps the integer power semantics of the FMM
  code, `freethreading_compatible=True`.
- The backward transform with `method=0` was never usable before this release; its first-order accuracy is
  low (relative error around 7e-2 on the Gaussian test case), which is expected for the method.

## Earlier versions

No changelog was kept before 0.7.0.
