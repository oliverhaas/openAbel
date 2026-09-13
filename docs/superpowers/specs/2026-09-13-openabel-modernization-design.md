# openAbel modernization: design

Date: 2026-09-13. Status: approved (sections 1 and the questions below in chat; the rest by "just do it").

## Goal

Replace openAbel's 2016-era tooling (`setup.py`-only build, Cython 0.29, nose, Travis, Sphinx on a
frozen Read the Docs site) with the tool set used in the author's other packages (uv, `pyproject.toml`,
Cython 3, pytest, ruff, ty, pre-commit, GitHub Actions, cibuildwheel, mkdocs-material + mike), fix the
bugs that make the current code unusable on a current Linux, and leave a PyPI release pipeline in place
but switched off. The numerical results of every working code path stay bit-for-bit or within
floating-point noise of today's.

## Decisions taken in the brainstorm

| Question | Decision |
|---|---|
| Publish now? | No. Repo stays private; the publish, tag and docs workflows are committed with a `.disabled` suffix. Activation checklist in section 6. |
| Code scope | Cython 3 migration, the allocation bug, the side-path bugs of section 2.3 (each with a regression test) and the even-order buffer bug found while writing those tests. No API change: parameter and function names stay camelCase (physics/maths code, not web code); ruff's naming rules are switched off for this repo. |
| Python floor / wheels | `>=3.12`; wheels for CPython 3.12, 3.13, 3.14 and 3.14t on manylinux x86_64 and macOS arm64. No Windows (`aligned_alloc` is missing from MSVC's CRT; fixing it touches every `free` site). |
| Docs | Port to mkdocs-material + mike, math through `pymdownx.arithmatex` + MathJax. `mkdocs build --strict` in CI. Site deployment dormant. |
| Approach | Template transplant from django-massless / django-filthyfields (the author's Cython repos), src layout. |
| Defaults | Version `0.7.0` (1.0 reserved for the public release); default branch renamed `master` -> `main` (done 2026-09-13); changelog at `docs/reference/changelog.md`; spec and plan under `docs/superpowers/`. |

Out of scope: openConv and openChargeState (each gets its own run later), any API redesign, lazy loading
of the coefficient tables, removing the `mathFun`/`complex.h` remnants, Windows support, Cython line-trace
coverage.

## 1. Repository layout and packaging

```
openAbel/
├── pyproject.toml  setup.py  MANIFEST.in  uv.lock  .gitignore  .pre-commit-config.yaml  mkdocs.yml
├── README.md (converted from README.rst)   LICENSE (unchanged, GPL-3.0)
├── src/openAbel/
│   ├── __init__.py                       exports Abel, __all__, __version__
│   ├── constants.{pxd,pyx}  helper.{pxd,pyx}  mathFun.{pxd,pyx}
│   └── abel/  __init__.py  coeffs.py  base/hansenLaw/trap/fmm/wrap.{pxd,pyx}  coeffsData/*.npy (238 files, 3.5 MB)
├── tests/  analytic.py  test_input.py  test_methods.py  test_package.py
├── docs/   index.md  transform-types.md  transform-methods.md  remarks.md  javascripts/  examples/  reference/  superpowers/
├── examples/ (kept as they are, formatted by ruff)          add/ (kept, excluded from lint)
└── .github/  dependabot.yml
    └── workflows/  ci.yml  dependabot-automerge.yml  publish.yml.disabled  tag.yml.disabled  docs.yml.disabled
```

Deleted: `.travis.yml`, `docs/conf.py`, `docs/Makefile`, every `docs/**/*.rst`, `README.rst`, the old
`setup.py` body. `add/` (Mathematica notebooks and the h5-to-npy converter) and `.gitattributes` stay.

### pyproject.toml

`[project]`: `name = "openAbel"`, `version = "0.7.0"` (single source of truth), `description`,
`readme = "README.md"`, `license = "GPL-3.0-or-later"`, `license-files = ["LICENSE"]` (PEP 639, so no
license classifier), `requires-python = ">=3.12"`, `dependencies = ["numpy>=2.0", "scipy>=1.13"]`,
`authors = [{ name = "Oliver Haas", email = "ohaas@e1plus.de" }]`, keywords (abel transform, fast multipole
method, end corrections, cython), classifiers: Development Status 4 - Beta, Intended Audience
Science/Research, Operating System POSIX :: Linux and MacOS, Programming Language Cython / Python 3 :: Only /
3.12 / 3.13 / 3.14, Topic Scientific/Engineering :: Physics and :: Mathematics. `[project.urls]`: Homepage
and Repository on GitHub; Documentation is added when the site exists.

`[build-system] requires = ["setuptools>=77", "cython>=3.1", "scipy>=1.13"]`, `build-backend =
"setuptools.build_meta"`. scipy is a build requirement because `cimport scipy.linalg.cython_blas` resolves
the `.pxd` at cythonize time (the extension then fetches the BLAS function pointers from scipy at import).

`[dependency-groups]`: `dev` = cython, setuptools, pytest, ruff, ty, pre-commit (exact `==` pins, as in the
other repos; dependabot keeps them current); `docs` = mkdocs, mkdocs-material, mike (exact pins).
`[tool.uv] default-groups = ["dev", "docs"]`, so a plain `uv sync` installs both (CI and the hooks rely on it).

Setuptools: `[tool.setuptools.packages.find] where = ["src"]`, `namespaces = false` (so `coeffsData/`,
which has no `__init__.py`, is data of `openAbel.abel` and not a namespace package);
`[tool.setuptools.package-data]` `openAbel = ["*.pxd"]`, `"openAbel.abel" = ["*.pxd", "coeffsData/*.npy"]`;
`[tool.setuptools.exclude-package-data] "*" = ["*.c", "*.html"]`.

`[tool.cibuildwheel]`: `build = ["cp312-*", "cp313-*", "cp314-*", "cp314t-*"]` (cibuildwheel >= 3.4 builds
free-threaded targets without an enable flag; the publish workflow pins `pypa/cibuildwheel@v3.4.1`), `build-frontend = "build[uv]"`, `build-verbosity = 1`,
`test-requires = ["pytest"]`, `test-command = "python -m pytest {project}/tests"`;
`[tool.cibuildwheel.linux] archs = ["x86_64"]`; `[tool.cibuildwheel.macos] archs = ["arm64"]`.

No `py.typed` and no `Typing :: Typed` classifier: the package has no type annotations, so the claim would
be false.

### setup.py

Only what setuptools cannot read from `pyproject.toml`: discover `src/openAbel/**/*.pyx`, map each to its
dotted module name (`src/openAbel/abel/fmm.pyx` -> `openAbel.abel.fmm`), build `Extension`s with
`extra_compile_args=["-O3"]` (the old `-Wunused-but-set-variable -Wsign-compare` enabled warnings rather
than suppressing them and are dropped), and call `cythonize(extensions, compiler_directives=...)` with
`language_level = "3"`, `boundscheck = False`, `nonecheck = False`, `wraparound = False`,
`cdivision = True`, `binding = True`, `cpow = True`, `freethreading_compatible = True`. `cpow = True`
keeps the Cython 0.29 semantics of `2 ** int` in `fmm.pyx` (Cython 3 would otherwise return a double).
The logo print goes.

### Package init and metadata

`src/openAbel/__init__.py`: `from .abel import Abel`, `__all__ = ["Abel"]`,
`__version__ = importlib.metadata.version("openAbel")`. `src/openAbel/abel/__init__.py` keeps
`from .wrap import Abel` with a `# ty: ignore[unresolved-import]` (ty cannot see compiled modules).

`MANIFEST.in`: `include LICENSE README.md`, `recursive-include src/openAbel *.pyx *.pxd *.npy`,
`prune tests`. This is what makes the sdist installable again (the 2020 sdist shipped `.c` files while
`setup.py` referenced `.pyx`).

`.gitignore`: the usual Python set plus `src/**/*.c`, `src/**/*.so`, `src/**/*.html`, `build/`, `dist/`,
`wheelhouse/`, `site/`, `.venv/`, `*.egg-info/`, `.claude/worktrees/`, `.superpowers/`.

Dev loop: `uv sync` builds the extensions in place (editable install); after editing a `.pyx`:
`uv sync --reinstall-package openAbel`.

## 2. Code changes

Every change below is small and local. Everything not listed stays byte-identical apart from `ruff format`
on the `.py` files (ruff never touches `.pyx`/`.pxd`) and pre-commit's `end-of-file-fixer`, which trims the
blank lines at the end of six `.pyx` files (there is no `trailing-whitespace` hook, so nothing else moves).

### 2.1 Cython 3 migration (build blockers and deprecations)

1. `abel/base.pyx:59`: the `else` branch of the method dispatch raises inside a nested `with gil:` while
   the GIL is already held; Cython 3 rejects that. Drop the inner `with gil:`, keep the `raise`. This edit
   lands with the packaging skeleton, because nothing builds without it.
2. `cpow = True` compiler directive (see setup.py above).
3. 24 signatures spelled `... nogil except -1` / `nogil except NULL` (in `base.pxd/.pyx`, `trap.pxd/.pyx`,
   `fmm.pxd/.pyx`) become `... except -1 nogil` / `except NULL nogil`; Cython 3 warns that the old order
   "will be disallowed".
4. `hansenLaw.pxd/.pyx`: the four `cdef int ...(...) nogil` entry points (`plan_fat_hansenLawOrgLin`,
   `plan_fat_hansenLawLinear`, `execute_fat_hansenLawLinear`, `destroy_fat_hansenLawLinear`) get an
   explicit `except -1`. Under Cython 0.29 the missing clause silently swallowed the `NotImplementedError`
   raised for the modified forward transform (`forwardBackward=-2`) and returned a copy of the input; under
   Cython 3 exceptions propagate by default, and the explicit clause documents that this is intended.
5. `freethreading_compatible = True` (setup.py) so importing on 3.14t does not re-enable the GIL.

### 2.2 Allocation (`helper.pxd`)

`cdef: size_t stdAlgn = max(sizeof(void*), 64)` in a `.pxd` is never executed, so every
`nullCheckMalloc`/`nullCheckCalloc` call has always passed `alignment = 0` to `aligned_alloc`. glibc <=
2.37 silently treated that as `malloc`; glibc >= 2.38 (Ubuntu 24.04) returns NULL, and the "exit(-1)" that
follows is Python's `builtins.exit` raising `SystemExit` inside a `nogil` function, which is written as an
unraisable and ignored, so callers proceed to write through NULL. Result today: every `Abel()` construction
dies with exit status 255 on this machine.

New `helper.pxd`:

- `stdAlgn` is removed; both functions default to `size_t alignment = 64`.
- The requested size is rounded up to a multiple of `alignment` (C11 requires it; macOS's libmalloc
  enforces it), and a zero-byte request allocates one `alignment` block so the returned pointer is always
  valid and freeable with `free()`.
- On NULL the functions raise `MemoryError` (`with gil: raise MemoryError(...)`) instead of printing and
  calling `exit`. Both are declared `except NULL nogil` (the form `base.pyx` already uses for `plan_fat`);
  the callers' existing `try/except` cleanup paths then run. `nullCheckCalloc` keeps the `memset`.
- Names and the `cimport ... nullCheckMalloc as malloc, nullCheckCalloc as calloc` call sites stay.

### 2.3 Side-path bugs (regression tests in section 3)

| # | Where | Bug | Fix |
|---|---|---|---|
| a | `abel/trap.pyx:93` | Backward transform with `method=0` loads `.../data/coeffs_deriv_smooth_02.npy`; the directory is `coeffsData/` and the loader for it is `coeffs.getCoeffs`. Every `Abel(n, 1, ..., method=0)` raises `FileNotFoundError`. | `coeffs.getCoeffs("coeffs_deriv_smooth", 2)` (shape `(3,)`, matches `orderFilter = 3`). Drop the now-unused `os.path` and `datetime` imports. |
| b | `abel/trap.pyx:410`, `abel/fmm.pyx:180` | Key `"coeffs_invSqrtDiffSqY2OR2_sing_small_halfShift_"` (trailing underscore) never exists; `forwardBackward=-2, shift=0.5` raises `KeyError` with `method=2` and crashes with `method=3`. | Key `"coeffs_invSqrtDiffSqY2OR2_sing_small_halfShift"`. |
| c | `abel/fmm.pyx:58-59` | `md.ltp` and `md.direct0` are not in the NULL-initialisation list but are freed by `destroy_fat_fmmTrapEndCorr`; any exception raised while planning before they are allocated (unsupported `shift`, bug b, `MemoryError`) frees garbage: segfault or `free(): invalid pointer`. Reproduced with `Abel(200, -1, 0.25, 0.01, method=3)`. | Add `md.ltp = md.direct0 = NULL` to the list. |
| d | `abel/hansenLaw.pyx:53,109` | See 2.1 item 4. The raise paths also leak: `plan_fat_hansenLawLinear` leaves `plan.methodData` allocated and `execute_fat_hansenLawLinear` leaks `xk`. | `except -1` on the four functions, `NotImplementedError("Method not implemented for given parameters.")` messages, `destroy_fat_hansenLawLinear(plan)` before the raise in plan, `free(xk)` before the raise in execute. |
| e | `abel/trap.pyx` (`execute_fat_trapezoidalEndCorr`), `abel/fmm.pyx` (`execute_fat_fmmTrapEndCorr`) | Found while writing the tests. The temporary input buffers are sized `N+order+orderFilter-2` and `N+order-1`, but the stencil reach per side is `(order-1)//2` plus `(orderFilter-1)//2`: for even `order` (the default 2 included) the buffers are one element too long, the copy loop reads `dataIn[N]` past the end of the input, and the FMM's direct-block `dgemv` multiplies the extra (with boundary 0: uninitialised) element of `dataInTemp1` by a zero coefficient. A NaN or Inf bit pattern there poisons the result; the `fb=1, method=3` cell failed sporadically in the probe runs. | Buffers `N+2*(ordM1Hlf+ordFilM1Hlf)` and `N+2*ordM1Hlf`, copy count and convolve length to match, `dgemv` row count clamped to `min(N-kk*ss-1, 2*ss)` (the rows beyond the data end are zero in `md.direct`, so the results are unchanged). Checked with AddressSanitizer. |
| f | `abel/trap.pyx` (`plan_fat_trapezoidalDesingConst`) | `md` is allocated in the `cdef` block and again in the body; the first block leaks on every `method=0` construction. | Declare `md` without the initialiser. No test: a leak is not observable from pytest. |

### 2.4 `abel/coeffs.py`

Rewritten with `pathlib` and annotations, same behaviour: eager load of every `coeffsData/*.npy` at import
into `{family: {order: array}}`, `getCoeffs(coeffsName, order)` unchanged in name and semantics. Family =
file name without `_NN.npy`, order = the two-digit suffix.

## 3. Tests

pytest; `nose` goes. `tests/` is a plain directory (no `__init__.py`), `testpaths = ["tests"]`,
`xfail_strict = true`. No coverage plugin: the code under test is Cython and would need line-tracing builds
to be measured; a coverage number for the 40 lines of Python would mislead.

- `tests/analytic.py` (helper, not collected): the Gaussian `exp(-x^2)` test pair for every transform type
  on the grid `x_i = (i + shift) * stepSize`, with the truncated-domain references at `R = x[-1]`:
  forward `sqrt(pi) exp(-y^2) erf(sqrt(R^2 - y^2))`; backward and backward-with-derivative
  `exp(-r^2) erf(sqrt(R^2 - r^2))` from inputs `sqrt(pi) exp(-y^2)` and `-2 y sqrt(pi) exp(-y^2)`; modified
  forward `y^2 exp(-y^2) [k1e(y^2/2) - k0e(y^2/2)]` (2 at y = 0) minus the tail
  `2 y^2 int_R^inf exp(-r^2) / (r^2 sqrt(r^2 - y^2)) dr` by `scipy.integrate.quad`. These are the
  formulas of the 1152-combination baseline sweep. Keyword-only API: `grid(*, shift)`,
  `relativeError(*, dataOut, reference)`, `inputSamples(*, forwardBackward, x)` (the input function at
  arbitrary points, also outside `[0, R]`), `analyticPair(*, forwardBackward, shift)`; constants `N_DATA = 200`,
  `X_MAX = 3.5`, `STEP_SIZE = X_MAX / (N_DATA - 1)`.
- `tests/test_methods.py`: behaviour preservation. `N_DATA = 200`, `X_MAX = 3.5`, error metric
  `max |out - ref| / max |ref|` over all but the last sample (the last sample is 0 by construction and is
  asserted separately). Parametrised over `forwardBackward in (-1, 1, 2, -2)`, `shift in (0.0, 0.5)`,
  `method in (2, 3)`, `order in (1, 2, 3, 5, 10)` with the tolerance table below, plus `method=0` and
  `method=1` at their single order. Tolerances are the larger of the two shifts' errors measured on the
  Cython 3 probe build (2026-09-13), times 5, rounded up to the next power of ten, so they are at least 6x
  above the observed value and far above any BLAS/platform noise (~1e-15). Three kinds of cells deviate: the
  `o10` cells sit 10-100x above the rule (the measured errors are 1e-15 to 7e-13, within reach of BLAS noise),
  and `fb=2, method=1` and `fb=1, method=0` are `2e-1` (three times the measured `6.4e-2` and `6.8e-2`; the
  rule would give 1):

  | fb | o1 | o2 | o3 | o5 | o10 | m0 | m1 |
  |---|---|---|---|---|---|---|---|
  | -1 | 1e-2 | 1e-4 | 1e-6 | 1e-9 | 1e-13 | 1e-2 | 1e-2 |
  | 1 | 1e-1 | 1e-1 | 1e-4 | 1e-6 | 1e-11 | 2e-1 | 1e-1 |
  | 2 | 1e-1 | 1e-4 | 1e-4 | 1e-7 | 1e-13 | 1e-1 | 2e-1 |
  | -2 | 1e-2 | 1e-3 | 1e-6 | 1e-9 | 1e-11 | 1e-2 | raises NotImplementedError |

  (`fb=1, o2` is `1e-1` because the `shift=0.5` error is `1.14e-2`, larger than at `o1`; `shift=0` gives
  `2.2e-4`.) Two combinations never ran before the fixes: `fb=1, method=0` (bug a) measures `6.8e-2` at
  `shift=0` and `9.5e-3` at `shift=0.5` (first-order accuracy, expected for the method), hence `2e-1`;
  `fb=-2, shift=0.5` with `method=2` and `method=3` (bug b) satisfies the `fb=-2` row as it stands.
- Boundary value 3: the same end-correction grid once more with `leftBoundary=rightBoundary=3`, the input
  extended by `outsideSamplesPerSide(*, forwardBackward, order) = (order-1)//2 + (orderFilter-1)//2` samples per
  side (`orderFilter = order+1+order%2` for `fb=1`, else 1) taken from `analytic.inputSamples`, same
  tolerances. This is the path bug e's buffer bookkeeping changes most.
- Regression tests live in the module that owns the subject (no `test_regressions.py`), the bug context in a
  one-line comment: (a) the `fb=1, method=0` cells of the single-order grid in `test_methods.py` (output
  finite and within tolerance); (b) `test_modifiedForwardHalfShift_trapezoidalAndFmmAgree` in
  `test_methods.py`: `fb=-2, shift=0.5` with `method=2` and `method=3` agree to `1e-8` for
  `order in (1, 2, 3, 5, 10)`, and the grid cells cover the `fb=-2` tolerances; (c)
  `test_unsupportedShift_raisesNotImplemented` in `test_input.py`: `Abel(200, fb, 0.25, 0.01, method)` for
  every `fb` and `method in (2, 3)` raises `NotImplementedError` (no crash); (d)
  `test_hansenLawModifiedForward_raisesNotImplemented` and `test_hansenLaw_worksAfterFailedConstruction` in
  `test_methods.py`: `Abel(200, -2, shift, step, method=1)` raises `NotImplementedError` and
  `Abel(200, -1, 0.0, step, method=1).execute(...)` still works afterwards; (e)
  `test_fmmBackwardOrder2_ignoresInputBeyondStencilReach` in `test_methods.py`: a NaN planted one sample
  beyond what boundary value 3 consumes must not reach the output (fails before the fix, passes after).
- `tests/test_input.py`: `method in (-1, 5)` raise `NotImplementedError` (the old file defined the same
  test name twice, so `-1` was never tested); `order=0` raises `ValueError` with `method in (2, 3)`; unsupported
  `shift` with `method in (2, 3)` raises `NotImplementedError` (test c above).
- `tests/test_package.py`: `openAbel.__version__` matches `\d+\.\d+\.\d+`, `openAbel.__all__ == ["Abel"]`.

197 tests in total (2 package, 12 input, 183 methods); the suite runs on 3.12, 3.13, 3.14 and 3.14t locally,
and once under AddressSanitizer for bug e.

The cibuildwheel `test-command` runs this same suite against every built wheel.

## 4. Lint, types, hooks

ruff (`pyproject.toml`): `target-version = "py312"`, `line-length = 120`, `fix = true`, `select = ["ALL"]`,
`ignore = ["COM812", "CPY", "D", "E501", "EM", "N", "TRY003"]` (`N`: camelCase API and module name are
deliberate; `CPY`: no copyright headers), `extend-exclude = ["add"]`. Per-file: `tests/**` = `["ANN", "INP001",
"PLR2004", "PT011", "S101"]`; `examples/**` = `["ANN", "ARG001", "B007", "ERA001", "F841", "FBT003", "ICN001",
"INP001", "NPY002", "PLR2004", "PTH", "T201"]` (the examples are kept as written; the ignores cover what
`ruff --fix` cannot fix in them). ruff sees only `.py` files (`__init__.py`, `coeffs.py`, `setup.py`, tests,
examples); never pass `.pyx` paths to it. The hooks reformat the seven example scripts once.

ty: `uv run ty check src tests`; the single compiled-import site carries an inline ignore (section 1). No
other overrides are needed.

pre-commit: the django-massless configuration (pre-commit-hooks v5.0.0, taplo v0.9.3 with `--reorder-arrays
--reorder-keys`, add-trailing-comma v3.1.0, local `uv sync`, `ruff check --fix`, `ruff format`) with
`uv run ty check src tests` as the pre-push hook, `default_language_version: python3.12` and `exclude: ^add/`
(the Mathematica exports). `end-of-file-fixer` also trims the trailing blank lines of six `.pyx` files
(section 2). The Sphinx tree is removed before the hooks first run: `docs/conf.py` is a `.py` file ruff would
lint, and the hooks would touch the `.rst` files.

## 5. CI

`.github/workflows/ci.yml`, on push to `main` and pull requests to `main`:

- `lint` (ubuntu-latest, `UV_PYTHON: "3.14"`): `uv sync` (default groups `dev` and `docs`), `uv lock --check`,
  `ruff check --no-fix` (`pyproject.toml` sets `fix = true`, so a plain `ruff check` would fix and pass),
  `ruff format --check`, `ty check src tests`, `mkdocs build --strict`.
- `test`, `fail-fast: false`, matrix `os: ubuntu-latest` x `python: ["3.12", "3.13", "3.14", "3.14t"]`
  plus `os: macos-latest, python: "3.14"`; each leg `uv python install`, `uv sync` (this compiles the
  extensions) and `uv run pytest`, with `UV_PYTHON` set to the matrix version so the free-threaded leg cannot
  pick up the runner's system interpreter.
- Action pins: `actions/checkout@v7`, `astral-sh/setup-uv@v7` (dependabot bumps them).

`dependabot.yml`: `uv` and `github-actions`, weekly. `dependabot-automerge.yml`: `gh pr merge --auto
--squash` for dependabot PRs (as in the other repos).

## 6. Release pipeline (dormant)

Three workflows are committed with a `.disabled` suffix and the django-massless header explaining why a
rename, not a commented-out YAML, is used (GitHub parses commented-out `.yml` and shows failing runs):

- `publish.yml.disabled`: the django-filthyfields shape. On `v*` tags or `workflow_dispatch(version)`: a
  `test` job (ruff, ty, pytest), `build_wheels` on ubuntu-latest and macos-latest with
  `pypa/cibuildwheel@v3.4.1` (consumes `[tool.cibuildwheel]`), `build_sdist` with `uv build --sdist`, both
  uploading with `actions/upload-artifact@v7`, `publish` with `actions/download-artifact@v8`
  (`merge-multiple: true`) and `pypa/gh-action-pypi-publish@release/v1` under `environment: pypi` with
  `id-token: write` (Trusted Publishing, no token secret).
- `tag.yml.disabled`: after a green CI run on `main`, read the version with `uv version --short`, skip if the
  version is on PyPI or the tag exists, otherwise tag `v<version>` and `gh workflow run` for `publish.yml` and
  `docs.yml` (tags pushed with `GITHUB_TOKEN` do not trigger workflows on their own).
- `docs.yml.disabled`: the django-cachex shape; `mike deploy --push --update-aliases <version> latest` on
  tags, `mike deploy --push main` on `main`.

Activation checklist (also in the `publish.yml.disabled` header): make the repo public (or accept a private
release); register the Trusted Publisher on PyPI (project `openAbel`, owner `oliverhaas`, workflow
`publish.yml`, environment `pypi`); create the `pypi` environment in the GitHub repo settings; enable
GitHub Pages from the `gh-pages` branch; refresh the action pins in the three files (dependabot does not see
`.disabled` files); rename the three files; delete the Read the Docs project or leave a redirect; set
`[project.urls] Documentation`.

## 7. Documentation

`mkdocs.yml`: the django-cachex configuration (material theme with light/dark palettes, navigation
features, `pymdownx.highlight/inlinehilite/superfences`, `admonition`, `pymdownx.details`, `attr_list`,
`md_in_html`, `toc` with permalinks, `extra.version.provider: mike`, `exclude_docs: superpowers/`) plus
`pymdownx.arithmatex` with `generic: true`, `extra_javascript` = `javascripts/mathjax.js` (the standard
mkdocs-material MathJax loader) and the MathJax 3 CDN script, and `pymdownx.snippets` with
`base_path: [".", "docs"]` so pages can include `README.md` and `examples/*.py` from the repo root.

Pages (converted from the RST, content kept, headings fixed): `index.md` (includes `README.md` through a
snippet, so the README is the single source), `transform-types.md`, `transform-methods.md`, `remarks.md`,
`examples/index.md` and `examples/example000.md` ... `example005.md` (figure plus the script included by
snippet in a `{ .python }` fence, replacing `literalinclude`), `reference/api.md` (the `Abel.__init__` and
`Abel.execute` parameter documentation from `wrap.pyx`'s docstrings, hand-written because griffe cannot parse
`.pyx`, plus the number of outside samples boundary value 3 consumes), `reference/changelog.md` (`0.7.0`
entry listing sections 1-2 including bugs e and f; earlier versions: "no changelog was kept").
The six PNGs stay in `docs/examples/`. `nav` in that order.

`README.md`: converted from `README.rst`. Badges: CI only (the Travis and RTD badges go). Install section:
`pip install openAbel` "once released on PyPI"; until then `pip install git+https://github.com/oliverhaas/openAbel`;
development: `uv sync`, `uv run pytest`, `uv run pre-commit install`. Requirements: Python >= 3.12, Linux
or macOS. Remove the `sudo python setup.py install` and Ubuntu 16.04 remarks. Keep the physics description,
the method overview and the GPL notice. The README lands with the packaging skeleton: `pyproject.toml` names it
as `readme`, so the build needs it.

## 8. Sequencing and bookkeeping

1. `master` -> `main` rename: done.
2. Work happens on `feat/modern-tooling` in a worktree; spec and plan are committed there first.
3. Implementation order in the plan: packaging skeleton, src move and README (build must pass before
   anything else), Cython 3 signatures and the allocation fix, side-path bug fixes with the ported tests and
   their regression tests green on 3.12-3.14t, the buffer fix (bug e) with its test and an AddressSanitizer
   run, Sphinx removal plus lint/type/hooks, CI, dormant release workflows, docs port with the changelog,
   final full check (`uv build`, install of the sdist into a fresh venv, `python -m pytest`).
4. Commit messages: Conventional Commits, no attribution lines, no `git add -A`.
5. At the end: finishing-a-development-branch (merge into `main` is the user's call); in `~/hq`
   `phd/STATE.md` the "Toolchain builds again" item moves from Up next into the log with what was done, and
   the parked "Public release and licensing" note is corrected (openAbel has a GPL-3.0 `LICENSE`; the
   missing file is openChargeState's).

## Risks

- macOS is built only in CI (no macOS machine here). The `aligned_alloc` size rounding (2.2) is the one
  change made for it; if the macOS leg fails for another reason, the fix is made in a follow-up commit on
  the same branch, not by dropping the leg.
- Free-threaded wheels declare the modules import-safe without the GIL; nothing claims that one `Abel`
  object may be used from several threads at once, and the docs do not say so.
- Method 3 depends on scipy's BLAS at the 1e-15 level; the tolerances in section 3 are 6-60x above the
  measured errors, so BLAS differences cannot flip a test.
