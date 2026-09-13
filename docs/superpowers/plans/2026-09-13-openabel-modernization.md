# openAbel Modernization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace openAbel's 2016-era tooling (`setup.py`-only build, Cython 0.29, nose, Travis, Sphinx) with `pyproject.toml` + uv + a `src/` layout, Cython 3, pytest, ruff, ty, pre-commit, GitHub Actions and mkdocs-material at version 0.7.0, fixing the bugs that make the current code unusable on a current Linux, without changing the public API or the numerical results.

**Architecture:** The package moves to `src/openAbel/` unchanged apart from a handful of local Cython edits (exception-clause order, the allocator, five side-path bugs, one buffer-size bug); `setup.py` shrinks to the `cythonize` call and everything else lives in `pyproject.toml`. A pytest suite with analytic Gaussian references for all four transform types gates every code change. Tooling (pre-commit, CI, dormant release workflows, mkdocs-material docs) is transplanted from the author's other repositories.

**Tech Stack:** Python >=3.12, Cython 3.3.0, setuptools 84.0.0, numpy >=2.0, scipy >=1.13, uv, pytest 9.1.1, ruff 0.16.7, ty 0.0.80, pre-commit 4.6.2, mkdocs 1.6.1 + mkdocs-material 9.7.7 + mike 2.2.0, GitHub Actions, cibuildwheel 3.4.1.

**Spec:** `docs/superpowers/specs/2026-09-13-openabel-modernization-design.md`

## Global Constraints

- Python `>=3.12`. Wheels for CPython 3.12, 3.13, 3.14 and the free-threaded 3.14t on manylinux x86_64 and macOS arm64. No Windows.
- The public API does not change: `openAbel.Abel(nData, forwardBackward, shift, stepSize, method=3, order=2, eps=...)` and `Abel.execute(dataIn, leftBoundary=0, rightBoundary=0)`. Parameter, function and module names stay camelCase (physics/maths code); ruff's `N` rules are switched off and stay off. No API redesign, no renames.
- The numerical results of every working code path stay bit-for-bit or within floating-point noise of today's. The tolerance tables in `tests/test_methods.py` are the gate; never loosen a tolerance to make a test pass.
- Version is `0.7.0`, single source of truth `pyproject.toml` (`openAbel.__version__` reads it through `importlib.metadata`).
- Exact tool pins, do not "update" them: `cython==3.3.0`, `setuptools==84.0.0`, `pytest==9.1.1`, `ruff==0.16.7`, `ty==0.0.80`, `pre-commit==4.6.2`, `mkdocs==1.6.1`, `mkdocs-material==9.7.7`, `mike==2.2.0`; pre-commit hook revs `pre-commit-hooks v5.0.0`, `taplo-pre-commit v0.9.3`, `add-trailing-comma v3.1.0`; GitHub Actions `actions/checkout@v7`, `astral-sh/setup-uv@v7`, `actions/upload-artifact@v7`, `actions/download-artifact@v8`, `pypa/cibuildwheel@v3.4.1`, `pypa/gh-action-pypi-publish@release/v1`.
- The repository stays private. The publish, tag and docs workflows are committed with a `.disabled` suffix and are never renamed to `.yml` by this plan.
- The `.pyx`/`.pxd` edits are exactly the ones listed in Tasks 1-3, applied by the patch scripts given there; ruff and ty never see `.pyx`/`.pxd` files. Trailing whitespace inside the `.pyx` files stays (no `trailing-whitespace` hook); only end-of-file blank lines are trimmed (Task 4).
- `add/` (Mathematica exports) and `src/openAbel/abel/coeffsData/*.npy` (238 files) are never edited; the `.npy` count stays 238 in the repo, the sdist and the wheel.
- Commit messages follow Conventional Commits (`build:`, `fix:`, `chore:`, `style:`, `ci:`, `docs:`), one line, no body needed, **no attribution lines or trailers of any kind**. Stage files by name (`git add <file> ...`); never `git add -A`, `git add .`, `git commit -a`. Never `git stash`.
- uv only: `uv sync`, `uv run <cmd>`, `uv build`. Never bare `pip install`. After editing a `.pyx` or `.pxd`: `uv sync --reinstall-package openAbel`.
- Conventions for new Python code (tests and helpers): function-based pytest tests named `test_<what>_<expectedOutcome>`; regression tests live in the module that owns the subject (`tests/test_methods.py`, `tests/test_input.py`), with the bug context in a one-line `# Regression: ...` comment, never in a `test_regression*.py` file; keyword-only arguments (`def f(*, a: int) -> int`) and full type hints on new helper functions; one-line docstrings or none; module-level containers immutable (`MappingProxyType`, `tuple`); top-level imports; no `from __future__ import annotations`.

## Working method (applies to every task)

- Work in the worktree `/home/ohaas/phd/openAbel/.claude/worktrees/feat+modern-tooling` on branch `feat/modern-tooling`; run every command from that directory.
- Whole new files: write them with your file-writing tool, byte-for-byte as given in the task (the content blocks are verbatim; keep the single trailing newline).
- Cython edits: each task gives a small Python patch script. Write it to `.superpowers/patches/<name>.py` (create the directory; `.superpowers/` is git-ignored from Task 1 on and is never committed) and run it with `python3 .superpowers/patches/<name>.py`. Every script asserts the exact occurrence count of the text it replaces, so a failed assertion means the file is not in the expected state: stop and report, do not improvise an edit by hand.
- The scripts locate the repository root through `Path(__file__).resolve().parents[2]`, so they only work from that directory.
- `uv sync` output: compiler warnings about `-Wsign-compare` in the generated C are expected; a Cython or C *error* is not.
- Expected outputs given as `N passed in ...` are from `uv run pytest -q`; timings vary.

## File map (final tree; the task that creates or replaces each file)

```
openAbel/
├── pyproject.toml  setup.py  MANIFEST.in  uv.lock  .gitignore  README.md         (1)
├── .pre-commit-config.yaml                                                       (4)
├── mkdocs.yml                                                                    (7)
├── LICENSE  .gitattributes  add/                                                 (unchanged)
├── .github/dependabot.yml  workflows/ci.yml  workflows/dependabot-automerge.yml  (5)
├── .github/workflows/{publish,tag,docs}.yml.disabled                             (6)
├── src/openAbel/__init__.py  abel/__init__.py  abel/coeffs.py                    (1)
├── src/openAbel/helper.pxd                                                       (2)
├── src/openAbel/abel/{base,trap,fmm,hansenLaw}.{pxd,pyx}                         (1, 2, 3, 4: edits only)
├── src/openAbel/{constants,mathFun}.{pxd,pyx}  helper.pyx  abel/wrap.{pxd,pyx}   (moved in 1; EOF trim in 4)
├── src/openAbel/abel/coeffsData/*.npy (238)                                      (moved in 1)
├── tests/test_package.py                                                         (1)
├── tests/analytic.py  test_input.py  test_methods.py                             (2; test_methods.py extended in 3)
├── examples/*.py (7, reformatted)  examples/*.png (2)                            (4)
└── docs/index.md  transform-types.md  transform-methods.md  remarks.md
    docs/javascripts/mathjax.js  examples/index.md  examples/example000-005.md
    docs/reference/api.md  reference/changelog.md                                 (7)
    docs/examples/*.png (6)                                                       (unchanged)
    docs/superpowers/                                                             (spec + this plan)
```

Deleted: `README.rst` (1), `docs/Makefile`, `docs/conf.py`, `docs/*.rst`, `docs/examples/*.rst` (4), `.travis.yml` (5).

---

### Task 1: Packaging skeleton, src layout and the build-blocking Cython edit

**Files:**
- Create: `.gitignore`, `pyproject.toml`, `MANIFEST.in`, `README.md`, `tests/test_package.py`, `uv.lock` (generated)
- Replace: `setup.py`, `src/openAbel/__init__.py`, `src/openAbel/abel/__init__.py`, `src/openAbel/abel/coeffs.py`
- Move: `openAbel/` -> `src/openAbel/` (257 tracked files: 19 sources + 238 `.npy`)
- Modify: `src/openAbel/abel/base.pyx` (lines 57-61 of the original file)
- Delete: `README.rst`

**Interfaces:**
- Consumes: nothing (first task).
- Produces: an importable `openAbel` package built from `src/` with `openAbel.Abel`, `openAbel.__version__ == "0.7.0"`, `openAbel.__all__ == ["Abel"]`; `openAbel.abel.coeffs.getCoeffs(coeffsName: str, order: int) -> np.ndarray` (the `.pyx` modules call it as `coeffs.getCoeffs(...)` / `cffs.getCoeffs(...)`); the dev loop `uv sync`, `uv sync --reinstall-package openAbel`, `uv run pytest`; a `.gitignore` that covers `.superpowers/`, `src/**/*.c`, `src/**/*.so`, `dist/`, `site/`.
- Known state at the end of this task: `import openAbel` works, but constructing `openAbel.Abel(...)` still terminates the process with exit status 255 (the allocator bug, fixed in Task 2). The old nose tests in `tests/test_input.py` and `tests/test_methods.py` are not runnable (they import `nose`) and are replaced in Task 2; run only `tests/test_package.py` here.

- [ ] **Step 1: Move the package into `src/`**

```bash
mkdir -p src
git mv openAbel src/openAbel
git status --short | grep -c '^R'
```

Expected: `257`. `git status --short` shows lines such as `R  openAbel/abel/fmm.pyx -> src/openAbel/abel/fmm.pyx`.

- [ ] **Step 2: Write `.gitignore`**

```gitignore
# Python
__pycache__/
*.py[cod]
*.egg-info/
.venv/
.pytest_cache/
.ruff_cache/

# Cython build products (the .pyx/.pxd sources are what is tracked)
src/**/*.c
src/**/*.so
src/**/*.html
build/
dist/
wheelhouse/

# Docs
site/

# Local tooling state
.claude/worktrees/
.superpowers/
```

- [ ] **Step 3: Write `pyproject.toml`**

```toml
[build-system]
build-backend = "setuptools.build_meta"
requires = ["cython>=3.1", "scipy>=1.13", "setuptools>=77"]

[project]
authors = [{ email = "ohaas@e1plus.de", name = "Oliver Haas" }]
classifiers = [
  "Development Status :: 4 - Beta",
  "Intended Audience :: Science/Research",
  "Operating System :: MacOS",
  "Operating System :: POSIX :: Linux",
  "Programming Language :: Cython",
  "Programming Language :: Python :: 3 :: Only",
  "Programming Language :: Python :: 3.12",
  "Programming Language :: Python :: 3.13",
  "Programming Language :: Python :: 3.14",
  "Topic :: Scientific/Engineering :: Mathematics",
  "Topic :: Scientific/Engineering :: Physics",
]
dependencies = ["numpy>=2.0", "scipy>=1.13"]
description = "Fast Abel transforms of equispaced data: Fast Multipole Method with arbitrary-order end corrections, in Cython"
keywords = [
  "abel transform",
  "cython",
  "end corrections",
  "fast multipole method",
]
license = "GPL-3.0-or-later"
license-files = ["LICENSE"]
name = "openAbel"
readme = "README.md"
requires-python = ">=3.12"
version = "0.7.0"

[project.urls]
Homepage = "https://github.com/oliverhaas/openAbel"
Repository = "https://github.com/oliverhaas/openAbel"

[dependency-groups]
dev = [
  "cython==3.3.0",
  "pre-commit==4.6.2",
  "pytest==9.1.1",
  "ruff==0.16.7",
  "setuptools==84.0.0",
  "ty==0.0.80",
]
docs = ["mike==2.2.0", "mkdocs==1.6.1", "mkdocs-material==9.7.7"]

[tool.uv]
default-groups = ["dev", "docs"]

[tool.setuptools.packages.find]
namespaces = false
where = ["src"]

[tool.setuptools.package-data]
"openAbel" = ["*.pxd"]
"openAbel.abel" = ["*.pxd", "coeffsData/*.npy"]

[tool.setuptools.exclude-package-data]
"*" = ["*.c", "*.html"]

[tool.cibuildwheel]
build = ["cp312-*", "cp313-*", "cp314-*", "cp314t-*"]
build-frontend = "build[uv]"
build-verbosity = 1
test-command = "python -m pytest {project}/tests"
test-requires = ["pytest"]

[tool.cibuildwheel.linux]
archs = ["x86_64"]

[tool.cibuildwheel.macos]
archs = ["arm64"]

[tool.pytest.ini_options]
testpaths = ["tests"]
xfail_strict = true

[tool.ruff]
extend-exclude = ["add"]
fix = true
line-length = 120
target-version = "py312"

[tool.ruff.lint]
ignore = ["COM812", "CPY", "D", "E501", "EM", "N", "TRY003"]
select = ["ALL"]

[tool.ruff.lint.per-file-ignores]
"examples/**" = [
  "ANN",
  "ARG001",
  "B007",
  "ERA001",
  "F841",
  "FBT003",
  "ICN001",
  "INP001",
  "NPY002",
  "PLR2004",
  "PTH",
  "T201",
]
"tests/**" = ["ANN", "INP001", "PLR2004", "PT011", "S101"]

[tool.ty.environment]
python-version = "3.12"
```

- [ ] **Step 4: Replace `setup.py`**

Overwrite the whole file (the old one hand-lists extensions, prints a logo and says version 0.6):

```python
"""Build the Cython extension modules; everything else lives in pyproject.toml."""

from pathlib import Path

from Cython.Build import cythonize
from setuptools import Extension, setup

SRC = Path("src")
COMPILER_DIRECTIVES = {
    "binding": True,
    "boundscheck": False,
    "cdivision": True,
    "cpow": True,  # keep the Cython 0.29 semantics of `2 ** int` (integer result) in fmm.pyx
    "freethreading_compatible": True,
    "language_level": "3",
    "nonecheck": False,
    "wraparound": False,
}

extensions = [
    Extension(
        pyx.relative_to(SRC).with_suffix("").as_posix().replace("/", "."),
        [str(pyx)],
        extra_compile_args=["-O3"],
    )
    for pyx in sorted((SRC / "openAbel").rglob("*.pyx"))
]

setup(ext_modules=cythonize(extensions, compiler_directives=COMPILER_DIRECTIVES))
```

- [ ] **Step 5: Write `MANIFEST.in`**

```
include LICENSE README.md
recursive-include src/openAbel *.pyx *.pxd *.npy
prune tests
```

- [ ] **Step 6: Replace the package init modules and the coefficient loader**

`src/openAbel/__init__.py` (the old file is the single line `from .abel import Abel`):

```python
"""openAbel: fast Abel transforms of equispaced data."""

from importlib.metadata import version

from .abel import Abel

__all__ = ["Abel"]
__version__ = version("openAbel")
```

`src/openAbel/abel/__init__.py` (the old file is the single line `from .wrap import Abel`; ty cannot see the compiled module, hence the ignore):

```python
from .wrap import Abel  # ty: ignore[unresolved-import]

__all__ = ["Abel"]
```

`src/openAbel/abel/coeffs.py` (same behaviour as the old `os.listdir` loop: eager load of every `coeffsData/*.npy` into `{family: {order: array}}`, family = file name without the `_NN` suffix, order = the two-digit suffix; `getCoeffs` keeps its name and semantics):

```python
"""Precomputed end-correction and filter coefficients, loaded eagerly from ``coeffsData/*.npy`` at import."""

from pathlib import Path

import numpy as np

dataDir = Path(__file__).parent / "coeffsData"

# Outer key: coefficient family, i.e. the file name without its "_NN" suffix. Inner key: the order NN.
coeffsAllDict: dict[str, dict[int, np.ndarray]] = {}
for path in sorted(dataDir.glob("*.npy")):
    coeffsName, _, order = path.stem.rpartition("_")
    coeffsAllDict.setdefault(coeffsName, {})[int(order)] = np.load(path).astype(np.double)


def getCoeffs(coeffsName: str, order: int) -> np.ndarray:
    """Return the coefficients of family ``coeffsName`` for the given ``order``."""
    return coeffsAllDict[coeffsName][order]
```

- [ ] **Step 7: Apply the one Cython 3 build blocker in `base.pyx`**

`plan_fat` raises inside a nested `with gil:` while the GIL is already held (the whole `try` sits in a `with gil:` block); Cython 3 rejects that. Write `.superpowers/patches/task1_base_with_gil.py`:

```python
"""Drop the nested `with gil:` around the raise in plan_fat's else branch (Cython 3 rejects it)."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
path = ROOT / "src/openAbel/abel/base.pyx"
source = path.read_text()
old = """            else:
                with gil:
                    raise NotImplementedError('Method not implemented for given parameters.')
        except:"""
new = """            else:
                raise NotImplementedError('Method not implemented for given parameters.')
        except:"""
assert source.count(old) == 1
path.write_text(source.replace(old, new))
print("patched", path.relative_to(ROOT))
```

Run: `mkdir -p .superpowers/patches && python3 .superpowers/patches/task1_base_with_gil.py`
Expected: `patched src/openAbel/abel/base.pyx`; `git diff --stat src/` reports 1 file changed, 1 insertion, 2 deletions.

- [ ] **Step 8: Replace `README.rst` with `README.md`**

```bash
git rm -q README.rst
```

Write `README.md`:

````markdown
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
pip install openAbel
```

Until then, install from the repository (this compiles the Cython extensions, so a C compiler is needed):

```bash
pip install git+https://github.com/oliverhaas/openAbel
```

A forward transform of a Gaussian sampled on 200 points, the first sample at `x = 0`:

```python
import numpy as np
import openAbel

nData = 200
stepSize = 3.5 / (nData - 1)
x = np.arange(nData) * stepSize

abelObj = openAbel.Abel(nData, -1, 0.0, stepSize)
dataOut = abelObj.execute(np.exp(-(x**2)))
```

The [examples](https://github.com/oliverhaas/openAbel/tree/main/examples) show the transform types and
methods in more detail, starting with `example000_simpleForward.py`; they need matplotlib.

## Development

```bash
uv sync                    # builds the extensions into .venv (editable install)
uv run pytest
uv run pre-commit install  # ruff and the other hooks on commit, ty on push
```

After editing a `.pyx` or `.pxd` file, rebuild with `uv sync --reinstall-package openAbel`.

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
coefficients up to 20th order, otherwise it's recommended to use at most 5th order. The FMM leads to a
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
````

(The CI badge points at the workflow Task 5 adds; the docs links point at the pages Task 7 adds.)

- [ ] **Step 9: Build and import**

```bash
uv sync --python 3.14
uv run python -c "import openAbel; print(openAbel.__version__, openAbel.__all__)"
```

Expected: `uv sync` creates `.venv` with CPython 3.14, resolves the dependencies, writes `uv.lock`, builds `openabel` from the worktree (two to three minutes) and ends with the installed-package list containing `+ openabel==0.7.0 (from file:///...)`. The python one-liner prints `0.7.0 ['Abel']`.

- [ ] **Step 10: Write `tests/test_package.py` and run it**

```python
import re

import openAbel


def test_version_isSemver():
    assert re.fullmatch(r"\d+\.\d+\.\d+", openAbel.__version__)


def test_publicApi_isAbelOnly():
    assert openAbel.__all__ == ["Abel"]
```

Run: `uv run pytest tests/test_package.py -q`
Expected: `2 passed in ...`

- [ ] **Step 11: Confirm the known failure state (do not fix it here)**

```bash
uv run python -c "import openAbel; openAbel.Abel(200, -1, 0.0, 0.01)"; echo "exit status $?"
```

Expected: `exit status 255` (an allocation error message may precede it). This is the `aligned_alloc(0, ...)` bug that Task 2 fixes; leave it.

- [ ] **Step 12: Commit**

`git status --short` must list only: the 257 renames, `M src/openAbel/__init__.py`, `M src/openAbel/abel/__init__.py`, `M src/openAbel/abel/coeffs.py`, `M src/openAbel/abel/base.pyx`, `M setup.py`, `D README.rst`, and the untracked `.gitignore`, `pyproject.toml`, `MANIFEST.in`, `README.md`, `uv.lock`, `tests/test_package.py`. (`.superpowers/`, `.venv/`, `src/**/*.c`, `src/**/*.so` and `*.egg-info/` are ignored and must not appear.)

```bash
git add .gitignore pyproject.toml setup.py MANIFEST.in README.md uv.lock tests/test_package.py
git add src/openAbel/__init__.py src/openAbel/abel/__init__.py src/openAbel/abel/coeffs.py src/openAbel/abel/base.pyx
git commit -m "build: switch to pyproject.toml, uv and a src layout"
```

(The renames were staged by `git mv`.) `git status --short` is empty afterwards.

---

### Task 2: Cython 3 signatures, the allocator, the side-path bugs and the pytest suite

Three commits, in this order. The full test grid cannot run before all three land: before the allocator fix every construction exits 255, and before bug (c) is fixed one grid cell aborts the whole pytest process.

**Files:**
- Modify: `src/openAbel/abel/base.pxd`, `base.pyx`, `trap.pxd`, `trap.pyx`, `fmm.pxd`, `fmm.pyx`, `hansenLaw.pxd`, `hansenLaw.pyx` (patch scripts)
- Replace: `src/openAbel/helper.pxd`
- Create: `tests/analytic.py`
- Replace: `tests/test_input.py`, `tests/test_methods.py` (the old nose files)

**Interfaces:**
- Consumes: the built package and `coeffs.getCoeffs` from Task 1.
- Produces: C-level `nullCheckMalloc(size_t MemSize, size_t alignment=64) except NULL nogil` and `nullCheckCalloc(size_t nn, size_t size, size_t alignment=64) except NULL nogil` in `helper.pxd` (call sites keep `cimport ... nullCheckMalloc as malloc, nullCheckCalloc as calloc`); every C-level entry point in `base/trap/fmm/hansenLaw` declared `except -1 nogil` (or `except NULL nogil`). Test support module `tests/analytic.py` with `N_DATA = 200`, `X_MAX = 3.5`, `STEP_SIZE`, `grid(*, shift: float) -> np.ndarray`, `relativeError(*, dataOut, reference) -> float`, `inputSamples(*, forwardBackward: int, x: np.ndarray) -> np.ndarray`, `analyticPair(*, forwardBackward: int, shift: float) -> tuple[np.ndarray, np.ndarray]`; `tests/test_methods.py` with the immutable tolerance tables `END_CORRECTION_TOLERANCE[(forwardBackward, order)]` and `SINGLE_ORDER_TOLERANCE[(forwardBackward, method)]`. Task 3 extends `tests/test_methods.py`.

#### Commit 1 of 3: except clauses after `nogil`

- [ ] **Step 1: Write and run the reorder script**

Cython 3 wants `except -1 nogil`; the old order `nogil except -1` still compiles but Cython warns it "will be disallowed". The four Hansen-Law entry points have no except clause at all; under Cython 3 exceptions propagate from them anyway, and the explicit clause documents that. Write `.superpowers/patches/task2_except_nogil.py`:

```python
"""Move the except clause after nogil in every C-level signature (Cython 3 order); make the Hansen-Law clauses explicit."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ABEL = ROOT / "src/openAbel/abel"

# (file, old, new, expected occurrence count)
EDITS = (
    ("base.pxd", "nogil except NULL", "except NULL nogil", 1),
    ("base.pyx", "nogil except NULL", "except NULL nogil", 1),
    ("base.pxd", "nogil except -1", "except -1 nogil", 2),
    ("base.pyx", "nogil except -1", "except -1 nogil", 2),
    ("trap.pxd", "nogil except -1", "except -1 nogil", 6),
    ("trap.pyx", "nogil except -1", "except -1 nogil", 6),
    ("fmm.pxd", "nogil except -1", "except -1 nogil", 3),
    ("fmm.pyx", "nogil except -1", "except -1 nogil", 3),
    ("hansenLaw.pxd", ") nogil\n", ") except -1 nogil\n", 3),
    ("hansenLaw.pyx", ") nogil:\n", ") except -1 nogil:\n", 4),
)
for name, old, new, count in EDITS:
    path = ABEL / name
    source = path.read_text()
    assert source.count(old) == count, (name, old, source.count(old))
    path.write_text(source.replace(old, new))
    print(f"{name}: {count} x {old.strip()!r} -> {new.strip()!r}")
```

Run: `python3 .superpowers/patches/task2_except_nogil.py`
Expected: ten lines, one per edit, no assertion error. Then:

```bash
grep -rn 'nogil except' src/openAbel; echo "remaining: $?"
grep -c 'except -1 nogil' src/openAbel/abel/hansenLaw.pxd src/openAbel/abel/hansenLaw.pyx
```

Expected: `remaining: 1` (grep found nothing) and `hansenLaw.pxd:3`, `hansenLaw.pyx:4`.

- [ ] **Step 2: Rebuild and import**

```bash
uv sync --reinstall-package openAbel
uv run python -c "import openAbel; print(openAbel.__version__)"
```

Expected: build succeeds, prints `0.7.0`.

- [ ] **Step 3: Commit**

```bash
git add src/openAbel/abel/base.pxd src/openAbel/abel/base.pyx src/openAbel/abel/trap.pxd src/openAbel/abel/trap.pyx src/openAbel/abel/fmm.pxd src/openAbel/abel/fmm.pyx src/openAbel/abel/hansenLaw.pxd src/openAbel/abel/hansenLaw.pyx
git commit -m "build: put except clauses after nogil for Cython 3"
```

#### Commit 2 of 3: the allocator

Background: the old `helper.pxd` declared `cdef: size_t stdAlgn = max(sizeof(void*), 64)` in the `.pxd`, which is never executed, so every allocation passed alignment 0 to `aligned_alloc`. glibc >= 2.38 returns NULL for that, and the `exit(-1)` that followed is Python's `builtins.exit` raising `SystemExit` inside a `nogil` function, written as an unraisable and ignored; callers then wrote through NULL. Hence exit status 255 on every construction.

- [ ] **Step 4: Replace `src/openAbel/helper.pxd`**

Overwrite the whole file with:

```cython
from libc.string cimport memset
cdef extern from "stdlib.h":
    void* aligned_alloc(size_t alignment, size_t size) nogil


# Inline malloc with null check, kind of a simple "hacky" solution not to have to do null check every time manually.
# Good enough for me and saves a lot of lines.
# Just "cimport [...] nullCheckMalloc as malloc" to replace normal malloc
# https://stackoverflow.com/questions/26831981/should-i-check-if-malloc-was-successful/26844703
# I decided to force alignment for up to AVX512 here, since it's usually worth it and not much lost if not.
# Might change this in the future. So for very specific cases alignment should be chosen manually anyway.
#
# C11 requires the size passed to aligned_alloc to be a multiple of the alignment (macOS enforces this, glibc does
# not), so the size is rounded up, and a zero-byte request allocates one alignment block so that the returned pointer
# is always valid and can be passed to free(). On failure a MemoryError is raised (the functions are "except NULL"),
# so the callers' cleanup paths run instead of writing through a NULL pointer.


cdef inline void* nullCheckMalloc(size_t MemSize, size_t alignment = 64) except NULL nogil:

    cdef:
        size_t allocSize = ((MemSize + alignment - 1) // alignment) * alignment
        void* AllocMem

    if allocSize == 0:
        allocSize = alignment

    AllocMem = aligned_alloc(alignment, allocSize)

    if NULL == AllocMem:
        with gil:
            raise MemoryError('aligned_alloc returned NULL, probably not enough memory or an invalid alignment.')

    return AllocMem


cdef inline void* nullCheckCalloc(size_t nn, size_t size, size_t alignment = 64) except NULL nogil:

    cdef:
        void* AllocMem = nullCheckMalloc(nn*size, alignment)

    memset(AllocMem, 0, nn*size)

    return AllocMem
```

- [ ] **Step 5: Rebuild and smoke-test a forward transform**

```bash
uv sync --reinstall-package openAbel
uv run python -c "import numpy as np, openAbel; x = np.arange(200) * 0.01; print(openAbel.Abel(200, -1, 0.0, 0.01).execute(np.exp(-x**2))[:3])"
```

Expected: three finite numbers close to `1.76` (the truncated forward transform of the Gaussian near `y = 0`), exit status 0. The 255 exit from Task 1 Step 11 is gone.

- [ ] **Step 6: Commit**

```bash
git add src/openAbel/helper.pxd
git commit -m "fix: request a valid alignment from aligned_alloc and raise MemoryError on NULL"
```

#### Commit 3 of 3: side-path bugs (a)-(d), the DesingConst double malloc, and the test suite

- [ ] **Step 7: Reproduce the bugs (red)**

Run each line separately and note the outcome.

```bash
uv run python -c "import openAbel; openAbel.Abel(200, 1, 0.0, 0.01, method=0)"
```
Expected (a): traceback ending in `FileNotFoundError` naming `.../data/coeffs_deriv_smooth_02.npy` (the directory is `coeffsData/`, and the loader for it is `coeffs.getCoeffs`).

```bash
uv run python -c "import openAbel; openAbel.Abel(200, -2, 0.5, 0.01, method=2)"
```
Expected (b): traceback ending in `KeyError: 'coeffs_invSqrtDiffSqY2OR2_sing_small_halfShift_'` (trailing underscore; the family is `coeffs_invSqrtDiffSqY2OR2_sing_small_halfShift`).

```bash
uv run python -c "import openAbel; openAbel.Abel(200, -1, 0.25, 0.01, method=3)"; echo "exit status $?"
```
Expected (c): the process aborts (`free(): invalid pointer`, `double free` or `Segmentation fault`; exit status 134 or 139) instead of ending in a `NotImplementedError` traceback: `md.ltp` and `md.direct0` are freed by the cleanup without ever having been NULL-initialised. (If it happens to end in the traceback, the garbage was NULL by chance; continue.)

```bash
uv run python -c "import openAbel; openAbel.Abel(200, -2, 0.0, 0.01, method=1)"
```
Expected (d): traceback ending in a bare `NotImplementedError` without message. Under Cython 0.29 this exception was swallowed and a copy of the input returned; under Cython 3 it propagates, but the plan data allocated before the raise leaks and the message is missing.

- [ ] **Step 8: Write and run the fix script**

Write `.superpowers/patches/task2_side_path_bugs.py`:

```python
"""Bugs a-d of spec section 2.3 plus the double malloc in plan_fat_trapezoidalDesingConst."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ABEL = ROOT / "src/openAbel/abel"

# (file, old, new, expected occurrence count)
EDITS = (
    # (a) backward transform with method 0: load the derivative filter through the coefficient loader
    (
        "trap.pyx",
        """                coeffs_filter_mv = np.load(osp.dirname(__file__) + "/data/coeffs_deriv_smooth_" + "%02d" % 2 + ".npy")\n""",
        """                coeffs_filter_mv = coeffs.getCoeffs('coeffs_deriv_smooth', 2)\n""",
        1,
    ),
    ("trap.pyx", "import os.path as osp\nimport datetime\n", "", 1),
    ("fmm.pyx", "import os.path as osp\nimport datetime\n", "", 1),
    # (b) the half-shift coefficient family has no trailing underscore
    (
        "trap.pyx",
        "'coeffs_invSqrtDiffSqY2OR2_sing_small_halfShift_'",
        "'coeffs_invSqrtDiffSqY2OR2_sing_small_halfShift'",
        1,
    ),
    (
        "fmm.pyx",
        "'coeffs_invSqrtDiffSqY2OR2_sing_small_halfShift_'",
        "'coeffs_invSqrtDiffSqY2OR2_sing_small_halfShift'",
        1,
    ),
    # (c) every pointer that destroy_fat_fmmTrapEndCorr frees starts out NULL
    (
        "fmm.pyx",
        "    md.direct = md.coeffsSing = md.coeffsNonsing = md.coeffsFilter = NULL\n",
        "    md.direct = md.direct0 = md.ltp = md.coeffsSing = md.coeffsNonsing = md.coeffsFilter = NULL\n",
        1,
    ),
    # (d) Hansen-Law: release what was allocated, then raise with a message
    (
        "hansenLaw.pyx",
        """    else:
        with gil:
            raise NotImplementedError

    return 0


# Hansen Law with linear approximation of function""",
        """    else:
        destroy_fat_hansenLawLinear(plan)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    return 0


# Hansen Law with linear approximation of function""",
        1,
    ),
    (
        "hansenLaw.pyx",
        """    else:
        with gil:
            raise NotImplementedError

    free(xk)""",
        """    else:
        free(xk)
        with gil:
            raise NotImplementedError('Method not implemented for given parameters.')

    free(xk)""",
        1,
    ),
    # plan_fat_trapezoidalDesingConst allocated md in the cdef block and again in the body, leaking the first block
    (
        "trap.pyx",
        "        methodData_DesingConst* md = <methodData_DesingConst*> malloc(sizeof(methodData_DesingConst))\n        int ii, jj, ll\n",
        "        methodData_DesingConst* md\n        int ii, jj, ll\n",
        1,
    ),
)
for name, old, new, count in EDITS:
    path = ABEL / name
    source = path.read_text()
    assert source.count(old) == count, (name, old, source.count(old))
    path.write_text(source.replace(old, new))
    print(f"{name}: {count} edit(s)")
```

Run: `python3 .superpowers/patches/task2_side_path_bugs.py`
Expected: nine lines (`trap.pyx: 1 edit(s)` x4, `fmm.pyx: 1 edit(s)` x3, `hansenLaw.pyx: 1 edit(s)` x2), no assertion error. `git diff --stat` reports 3 files changed.

Why the Hansen-Law cleanup is correct: `plan_fat` in `base.pyx` frees only `pl.grid` and `pl` on error, so the method-specific data must be released by the method itself; `destroy_fat_hansenLawLinear` frees `md.coeffs` and `md`. The DesingConst edit keeps the second `md = <methodData_DesingConst*> malloc(...)` in the function body (that is the one whose result is used).

- [ ] **Step 9: Rebuild and re-run the reproductions (green)**

```bash
uv sync --reinstall-package openAbel
uv run python -c "import openAbel; openAbel.Abel(200, 1, 0.0, 0.01, method=0); print('a ok')"
uv run python -c "import openAbel; openAbel.Abel(200, -2, 0.5, 0.01, method=2); print('b ok')"
uv run python -c "import openAbel; openAbel.Abel(200, -1, 0.25, 0.01, method=3)"; echo "exit status $?"
uv run python -c "import openAbel; openAbel.Abel(200, -2, 0.0, 0.01, method=1)"
```

Expected: `a ok`; `b ok`; a traceback ending in `NotImplementedError: Method not implemented for given parameters.` followed by `exit status 1` (no abort); a traceback ending in `NotImplementedError: Method not implemented for given parameters.`.

- [ ] **Step 10: Write `tests/analytic.py`**

```python
"""Analytic Gaussian test pairs for every openAbel transform type on the truncated domain [0, R]."""

import numpy as np
from scipy import integrate, special

# Test function f(x) = exp(-x^2) on the grid x_i = (i + shift) * STEP_SIZE. Every reference is the transform of the
# truncated integral with R = x[-1], which is what openAbel computes; the last sample is 0 by construction.
N_DATA = 200
X_MAX = 3.5
STEP_SIZE = X_MAX / (N_DATA - 1)


def grid(*, shift: float) -> np.ndarray:
    return (np.arange(N_DATA) + shift) * STEP_SIZE


def relativeError(*, dataOut: np.ndarray, reference: np.ndarray) -> float:
    """max |dataOut - reference| / max |reference| over all but the last sample."""
    return float(np.max(np.abs(dataOut[:-1] - reference[:-1])) / np.max(np.abs(reference)))


def _modifiedForwardTail(*, y: float, R: float) -> float:
    """2 y^2 int_R^inf exp(-r^2) / (r^2 sqrt(r^2 - y^2)) dr, the part of the modified forward transform beyond R."""
    if y == 0.0:
        return 0.0
    value, _ = integrate.quad(lambda r: np.exp(-(r**2)) / (r**2 * np.sqrt(r**2 - y**2)), R, np.inf, limit=200)
    return 2.0 * y**2 * value


def inputSamples(*, forwardBackward: int, x: np.ndarray) -> np.ndarray:
    """The input function of the transform type ``forwardBackward`` sampled at ``x`` (any points, also outside [0, R])."""
    g = np.exp(-(x**2))
    if forwardBackward in (-1, -2):
        # forward and modified forward transform f(r) = exp(-r^2)
        return g
    if forwardBackward == 1:
        # backward transform of F(y) = sqrt(pi) exp(-y^2); openAbel differentiates the input itself
        return np.sqrt(np.pi) * g
    if forwardBackward == 2:
        # backward transform with the derivative F'(y) = -2 y sqrt(pi) exp(-y^2) supplied as input
        return -2.0 * x * np.sqrt(np.pi) * g
    msg = f"No analytic input for forwardBackward={forwardBackward}"
    raise ValueError(msg)


def analyticPair(*, forwardBackward: int, shift: float) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(dataIn, expected dataOut)`` for the transform type ``forwardBackward`` on ``grid(shift)``."""
    x = grid(shift=shift)
    R = x[-1]
    g = np.exp(-(x**2))
    dataIn = inputSamples(forwardBackward=forwardBackward, x=x)
    truncatedErf = special.erf(np.sqrt(np.maximum(R**2 - x**2, 0.0)))
    if forwardBackward == -1:
        # forward: F(y) = 2 int_y^R r exp(-r^2) / sqrt(r^2 - y^2) dr = sqrt(pi) exp(-y^2) erf(sqrt(R^2 - y^2))
        return dataIn, np.sqrt(np.pi) * g * truncatedErf
    if forwardBackward in (1, 2):
        # backward: f(r) = exp(-r^2), truncated at R
        return dataIn, g * truncatedErf
    # modified forward: H(y) = 2 int_y^inf exp(-r^2) y^2 / (r^2 sqrt(r^2 - y^2)) dr
    #                        = y^2 exp(-y^2) [K1e(y^2/2) - K0e(y^2/2)], H(0) = 2, minus the tail beyond R
    h = np.full_like(x, 2.0)
    xp = x[x > 0.0]
    h[x > 0.0] = xp**2 * np.exp(-(xp**2)) * (special.k1e(0.5 * xp**2) - special.k0e(0.5 * xp**2))
    tail = np.array([*(_modifiedForwardTail(y=y, R=R) for y in x[:-1]), h[-1]])
    return dataIn, h - tail
```

- [ ] **Step 11: Replace `tests/test_input.py`**

The old file (nose) defined `test_methodNotImplemented` twice, so `method=-1` was never tested. Overwrite with:

```python
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
```

- [ ] **Step 12: Replace `tests/test_methods.py`**

Overwrite the old nose file with:

```python
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
```

- [ ] **Step 13: Run the suite**

Run: `uv run pytest -q`
Expected: `116 passed in ...` (2 package + 12 input + 102 methods). If exactly one `method=3` cell fails with a NaN or garbage result and passes when re-run, that is bug (e): an uninitialised buffer element multiplied by a zero coefficient (Task 3 fixes it); re-run once and continue. Any deterministic failure is a real problem: stop and report.

- [ ] **Step 14: Run the suite on every supported interpreter**

`uv sync --python X` recreates `.venv` for that interpreter (each run rebuilds the extensions, one to three minutes). Finish on 3.14 so that `.venv` ends on the interpreter CI lints with.

```bash
uv sync --python 3.12 --reinstall-package openAbel && uv run pytest -q
uv sync --python 3.13 --reinstall-package openAbel && uv run pytest -q
uv sync --python 3.14t --reinstall-package openAbel && uv run pytest -q
uv sync --python 3.14 --reinstall-package openAbel && uv run pytest -q
```

Expected: `116 passed` four times. `uv run python -c "import sys; print(sys.version)"` prints a 3.14 non-free-threaded version afterwards.

- [ ] **Step 15: Commit**

```bash
git add src/openAbel/abel/trap.pyx src/openAbel/abel/fmm.pyx src/openAbel/abel/hansenLaw.pyx tests/analytic.py tests/test_input.py tests/test_methods.py
git commit -m "fix: side-path bugs in the trapezoidal, FMM and Hansen-Law methods"
```

`git status --short` is empty afterwards (the extra `.so` files from the interpreter loop are ignored).

---

### Task 3: Exact end-correction input buffers (bug e)

**Files:**
- Modify: `src/openAbel/abel/trap.pyx` (`execute_fat_trapezoidalEndCorr`: buffer sizes, copy count, convolve length), `src/openAbel/abel/fmm.pyx` (`execute_fat_fmmTrapEndCorr`: the same three plus the direct-summation `dgemv` row count)
- Test: `tests/test_methods.py` (import line, one helper, two tests)

**Interfaces:**
- Consumes: `inputSamples(*, forwardBackward, x)`, `analyticPair`, `relativeError`, `N_DATA`, `STEP_SIZE` from `tests/analytic.py`; `END_CORRECTION_TOLERANCE` from `tests/test_methods.py` (Task 2).
- Produces: `outsideSamplesPerSide(*, forwardBackward: int, order: int) -> int` in `tests/test_methods.py`; the suite grows to 197 tests. Nothing later depends on new names.

Background: in both end-correction executes the temporary input buffers were sized `nData + order + orderFilter - 2` and `nData + order - 1`, while the stencil reach per side is `ordM1Hlf = (order-1)//2` plus `ordFilM1Hlf = (orderFilter-1)//2`. For odd `order` the two agree; for even `order` (including the default 2) the buffers are one element too long, the copy loop reads one sample past the end of `dataIn`, and the FMM's direct block `dgemv` feeds the extra (uninitialised, with boundary 0) element of `dataInTemp1` into the sum multiplied by a zero coefficient. A NaN or Inf bit pattern there poisons the result. The fix sizes the buffers `nData + 2*(ordM1Hlf + ordFilM1Hlf)` and `nData + 2*ordM1Hlf`, matches the copy count and the convolve length, and clamps the `dgemv` row count to the rows that hold data (`min(nData - kk*ss - 1, 2*ss)`; the rows beyond are zero in `md.direct`, so the results are unchanged).

- [ ] **Step 1: Write the failing regression test**

In `tests/test_methods.py`, change the import line

```python
from analytic import N_DATA, STEP_SIZE, analyticPair, relativeError
```

to

```python
from analytic import N_DATA, STEP_SIZE, analyticPair, inputSamples, relativeError
```

and insert, directly after `test_singleOrderMethods_matchAnalyticTransform` (before `test_modifiedForwardHalfShift_trapezoidalAndFmmAgree`):

```python
def outsideSamplesPerSide(*, forwardBackward: int, order: int) -> int:
    """Samples outside the domain that boundary value 3 consumes per side: the half widths of the end-correction
    stencil and, for the backward transform with numerical derivative, of the derivative filter."""
    orderFilter = order + 1 + order % 2 if forwardBackward == 1 else 1
    return (order - 1) // 2 + (orderFilter - 1) // 2


def test_fmmBackwardOrder2_ignoresInputBeyondStencilReach():
    # Regression: for even orders the end-correction methods read one input sample past the stencil reach, and the FMM
    # fed that sample (or, with boundary 0, an uninitialised buffer element) into its direct summation multiplied by a
    # zero coefficient. A NaN there poisoned the result; with boundary 0 it made the transform flaky.
    nOutside = outsideSamplesPerSide(forwardBackward=1, order=2)
    x = np.arange(-nOutside, N_DATA + nOutside + 1) * STEP_SIZE
    dataIn = inputSamples(forwardBackward=1, x=x)
    dataIn[-1] = np.nan  # one sample beyond what boundary value 3 needs
    abelObj = openAbel.Abel(N_DATA, 1, 0.0, STEP_SIZE, method=3, order=2)
    dataOut = abelObj.execute(dataIn, leftBoundary=3, rightBoundary=3)
    assert np.isfinite(dataOut).all()
```

- [ ] **Step 2: Run it to see it fail**

Run: `uv run pytest -q tests/test_methods.py -k ignoresInputBeyondStencilReach`
Expected: `1 failed` with `assert np.isfinite(dataOut).all()` -> `assert False` (the planted NaN reaches the output).

- [ ] **Step 3: Write and run the fix script**

Write `.superpowers/patches/task3_buffer_sizes.py`:

```python
"""Bug e: size the end-correction input buffers exactly and clamp the FMM direct dgemv to rows that hold data."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ABEL = ROOT / "src/openAbel/abel"

# (file, old, new, expected occurrence count)
EDITS = (
    (
        "trap.pyx",
        "    dataInTemp0 = <double*> malloc((pl.nData+md.order+md.orderFilter-2)*sizeof(double))\n    dataInTemp1 = <double*> malloc((pl.nData+md.order-1)*sizeof(double))\n",
        "    dataInTemp0 = <double*> malloc((pl.nData+2*(orderM1Half+orderFilterM1Half))*sizeof(double))\n    dataInTemp1 = <double*> malloc((pl.nData+2*orderM1Half)*sizeof(double))\n",
        1,
    ),
    (
        "trap.pyx",
        "    for ii in range(pl.nData+md.order+md.orderFilter-2-nLeftExt-nRightExt):\n",
        "    for ii in range(pl.nData+2*(orderM1Half+orderFilterM1Half)-nLeftExt-nRightExt):\n",
        1,
    ),
    (
        "trap.pyx",
        "    convolve(dataInTemp0, pl.nData+md.order-1, dataInTemp1, md.orderFilter, md.coeffsFilter)\n",
        "    convolve(dataInTemp0, pl.nData+2*orderM1Half, dataInTemp1, md.orderFilter, md.coeffsFilter)\n",
        1,
    ),
    (
        "fmm.pyx",
        "    dataInTemp0 = <double*> malloc((pl.nData+md.order+md.orderFilter-2)*sizeof(double))\n    dataInTemp1 = <double*> malloc((pl.nData+md.order-1)*sizeof(double))\n",
        "    dataInTemp0 = <double*> malloc((pl.nData+2*(ordM1Hlf+ordFilM1Hlf))*sizeof(double))\n    dataInTemp1 = <double*> malloc((pl.nData+2*ordM1Hlf)*sizeof(double))\n",
        1,
    ),
    (
        "fmm.pyx",
        "    for ii in range(pl.nData+md.order+md.orderFilter-2-nLeftExt-nRightExt):\n",
        "    for ii in range(pl.nData+2*(ordM1Hlf+ordFilM1Hlf)-nLeftExt-nRightExt):\n",
        1,
    ),
    (
        "fmm.pyx",
        "    convolve(dataInTemp0, pl.nData+md.order-1, dataInTemp1, md.orderFilter, md.coeffsFilter)\n",
        "    convolve(dataInTemp0, pl.nData+2*ordM1Hlf, dataInTemp1, md.orderFilter, md.coeffsFilter)\n",
        1,
    ),
    # the original dgemv line ends with "&mm, " (a trailing space before the newline)
    (
        "fmm.pyx",
        "    for kk in range(ll):\n        blas.dgemv('t', &mm, &md.ss, &ONED, &md.direct[kk*md.ss**2*2], &mm, \n",
        "    for kk in range(ll):\n        # Only the rows that hold data: rows beyond the data end are zero in md.direct, and reading the matching\n        # input elements would run past dataInTemp1.\n        nn = min(pl.nData-kk*md.ss-1, mm)\n        blas.dgemv('t', &nn, &md.ss, &ONED, &md.direct[kk*md.ss**2*2], &mm,\n",
        1,
    ),
)
for name, old, new, count in EDITS:
    path = ABEL / name
    source = path.read_text()
    assert source.count(old) == count, (name, old, source.count(old))
    path.write_text(source.replace(old, new))
    print(f"{name}: {count} edit(s)")
```

Run: `python3 .superpowers/patches/task3_buffer_sizes.py`
Expected: seven lines (`trap.pyx` x3, `fmm.pyx` x4), no assertion error. `git diff --stat src/` reports 2 files changed. The variables used already exist: `orderM1Half`/`orderFilterM1Half` are computed at the top of `execute_fat_trapezoidalEndCorr`, `ordM1Hlf`/`ordFilM1Hlf` and `nn` (declared `int`) at the top of `execute_fat_fmmTrapEndCorr`.

- [ ] **Step 4: Rebuild, run the regression test, then the whole suite**

```bash
uv sync --reinstall-package openAbel
uv run pytest -q tests/test_methods.py -k ignoresInputBeyondStencilReach
uv run pytest -q
```

Expected: `1 passed`, then `117 passed`.

- [ ] **Step 5: Add the boundary-3 grid test**

Append to `tests/test_methods.py`, directly after `test_fmmBackwardOrder2_ignoresInputBeyondStencilReach`:

```python
@pytest.mark.parametrize("method", [2, 3])
@pytest.mark.parametrize("shift", [0.0, 0.5])
@pytest.mark.parametrize(
    ("forwardBackward", "order", "tolerance"),
    [(forwardBackward, order, tolerance) for (forwardBackward, order), tolerance in END_CORRECTION_TOLERANCE.items()],
)
def test_endCorrectionMethods_outsideSamples_matchAnalyticTransform(forwardBackward, order, tolerance, shift, method):
    # Boundary value 3 on both sides: the input carries the samples the stencils reach into instead of extrapolating.
    nOutside = outsideSamplesPerSide(forwardBackward=forwardBackward, order=order)
    x = (np.arange(-nOutside, N_DATA + nOutside) + shift) * STEP_SIZE
    dataIn = inputSamples(forwardBackward=forwardBackward, x=x)
    _, reference = analyticPair(forwardBackward=forwardBackward, shift=shift)
    abelObj = openAbel.Abel(N_DATA, forwardBackward, shift, STEP_SIZE, method=method, order=order)
    dataOut = abelObj.execute(dataIn, leftBoundary=3, rightBoundary=3)
    assert dataOut.shape == (N_DATA,)
    assert dataOut[-1] == 0.0
    assert relativeError(dataOut=dataOut, reference=reference) < tolerance
```

Run: `uv run pytest -q`
Expected: `197 passed in ...`. Check the file shape: `grep -c '^def test_\|^def outsideSamplesPerSide' tests/test_methods.py` prints `8`.

- [ ] **Step 6: AddressSanitizer run**

Build the 3.14 extensions with ASan, run the suite with the sanitizer preloaded, then restore the normal build. `.venv/bin/python` must be called directly (a `uv run` wrapper would not carry `LD_PRELOAD` into the interpreter the same way).

```bash
CFLAGS="-fsanitize=address -fno-omit-frame-pointer" LDFLAGS="-fsanitize=address" uv sync --reinstall-package openAbel
nm -D src/openAbel/abel/fmm.cpython-314-x86_64-linux-gnu.so | grep -c asan
LD_PRELOAD=$(gcc -print-file-name=libasan.so) ASAN_OPTIONS=detect_leaks=0 .venv/bin/python -m pytest -q -p no:cacheprovider tests
uv sync --reinstall-package openAbel
nm -D src/openAbel/abel/fmm.cpython-314-x86_64-linux-gnu.so | grep -c asan
```

Expected: a positive count (about 30) after the ASan build; `197 passed` with no `ERROR: AddressSanitizer` block in the output; `0` after the restore. If `gcc -print-file-name=libasan.so` prints only `libasan.so` (library not installed) the preload fails to start: run `uv sync --reinstall-package openAbel` to restore the build, say so in the report and continue; the sanitizer run is a check, not a deliverable.

- [ ] **Step 7: Commit**

```bash
git add src/openAbel/abel/trap.pyx src/openAbel/abel/fmm.pyx tests/test_methods.py
git commit -m "fix: size the end-correction input buffers exactly"
```

---

### Task 4: Sphinx removal, pre-commit configuration and formatting

**Files:**
- Delete: `docs/Makefile`, `docs/conf.py`, `docs/examples.rst`, `docs/index.rst`, `docs/readmeLink.rst`, `docs/remarks.rst`, `docs/transformMethods.rst`, `docs/transformTypes.rst`, `docs/examples/example000.rst` ... `docs/examples/example005.rst`
- Create: `.pre-commit-config.yaml`
- Modify (by the hooks, no manual edits): `examples/example000_simpleForward.py`, `example001_simpleBackward.py`, `example002_methodOrder.py`, `example003_noisyBackward.py`, `example004_fullComparison.py`, `example005_comparisonPyAbel.py`, `example006_simpleForwardAndBackward.py` (ruff fixes, ruff format, trailing commas, end-of-file); `src/openAbel/abel/base.pyx`, `fmm.pyx`, `hansenLaw.pyx`, `trap.pyx`, `wrap.pyx`, `src/openAbel/mathFun.pyx` (end-of-file blank lines only)

**Interfaces:**
- Consumes: `pyproject.toml` `[tool.ruff]` and `[tool.ty]` (Task 1), the test suite (Tasks 2-3).
- Produces: the lint gate `uv run pre-commit run --all-files` that Tasks 5-8 run; `uv run ty check src tests` as the pre-push hook. Do **not** run `pre-commit install` during this plan (git hooks would then run on every later commit; the README tells developers to install them). `docs/` holds only the six PNGs and `docs/superpowers/` afterwards; Task 7 writes the mkdocs pages.

- [ ] **Step 1: Remove the Sphinx tree**

The Sphinx site is dead (the Read the Docs build has been frozen for years) and Task 7 replaces it. It has to go before the hooks run: `docs/conf.py` is a `.py` file that `ruff check` would lint with `select = ["ALL"]`, and `docs/Makefile` and five `.rst` files lack the single trailing newline `end-of-file-fixer` enforces.

```bash
git rm -q docs/Makefile docs/conf.py docs/examples.rst docs/index.rst docs/readmeLink.rst docs/remarks.rst docs/transformMethods.rst docs/transformTypes.rst
git rm -q docs/examples/example000.rst docs/examples/example001.rst docs/examples/example002.rst docs/examples/example003.rst docs/examples/example004.rst docs/examples/example005.rst
ls docs docs/examples
git commit -m "docs: remove the Sphinx configuration and pages"
```

Expected: `docs/` holds `examples/` and `superpowers/`; `docs/examples/` holds six PNGs (`example000_simpleForward.png` ... `example005_comparisonPyAbel.png`); `git status --short` is empty after the commit.

- [ ] **Step 2: Write `.pre-commit-config.yaml`**

```yaml
default_stages: [pre-commit]
default_install_hook_types:
  - pre-commit
  - pre-push
fail_fast: false
# Mathematica notebooks, kept as exported.
exclude: ^add/

default_language_version:
  python: python3.12

repos:
  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v5.0.0
    hooks:
      - id: check-ast
      - id: check-case-conflict
      - id: check-json
      - id: check-merge-conflict
      - id: check-symlinks
      - id: check-toml
      - id: check-yaml
      - id: debug-statements
      - id: detect-private-key
      - id: end-of-file-fixer
        stages: [pre-commit]
      - id: mixed-line-ending
        args: ["--fix=lf"]

  - repo: https://github.com/ComPWA/taplo-pre-commit
    rev: v0.9.3
    hooks:
      - id: taplo-format
        args: ["--", "--indent-string", "  ", "--reorder-arrays", "--reorder-keys"]

  - repo: https://github.com/asottile/add-trailing-comma
    rev: v3.1.0
    hooks:
      - id: add-trailing-comma

  - repo: local
    hooks:
      - id: uv-sync-check
        name: uv-sync-check
        language: system
        entry: uv sync
        pass_filenames: false

      - id: ruff-check
        name: ruff-check
        entry: uv run ruff check --fix
        language: system
        pass_filenames: false

      - id: ruff-format
        name: ruff-format
        entry: uv run ruff format
        language: system
        pass_filenames: false

  - repo: local
    hooks:
      - id: ty
        name: ty
        language: system
        entry: uv run ty check src tests
        pass_filenames: false
        always_run: true
        stages: [pre-push]
```

- [ ] **Step 3: Run the hooks until they are clean**

```bash
uv run pre-commit run --all-files
```

The first run installs the hook environments (network) and modifies files (`ruff-check`, `ruff-format`, `add-trailing-comma`, `end-of-file-fixer` report `Failed` with "files were modified by this hook"). Run the same command again, up to three times in total, until every hook reports `Passed` or `Skipped`. No manual edits and no `# noqa` are needed; if a hook keeps failing with a lint *error* (not a modification), stop and report. Note `ty` does not run here (pre-push stage).

- [ ] **Step 4: Check what changed**

```bash
git status --short
```

Expected: exactly these 13 modified files plus the untracked config:

```
 M examples/example000_simpleForward.py
 M examples/example001_simpleBackward.py
 M examples/example002_methodOrder.py
 M examples/example003_noisyBackward.py
 M examples/example004_fullComparison.py
 M examples/example005_comparisonPyAbel.py
 M examples/example006_simpleForwardAndBackward.py
 M src/openAbel/abel/base.pyx
 M src/openAbel/abel/fmm.pyx
 M src/openAbel/abel/hansenLaw.pyx
 M src/openAbel/abel/trap.pyx
 M src/openAbel/abel/wrap.pyx
 M src/openAbel/mathFun.pyx
?? .pre-commit-config.yaml
```

The Cython diffs must be removed blank lines at the end of each file and nothing else:

```bash
git diff src/ | grep '^[-+]' | grep -v '^[-+][-+]' | grep -v '^-$' | wc -l
```

Expected: `0`. (`pyproject.toml` is already in taplo's key order; `uv.lock` is unchanged by `uv sync`; `add/` is excluded.)

- [ ] **Step 5: Lint, type-check and test the result**

```bash
uv run ruff check --no-fix
uv run ruff format --check
uv run ty check src tests
uv sync --reinstall-package openAbel && uv run pytest -q
MPLBACKEND=Agg uv run --with matplotlib python examples/example000_simpleForward.py && rm -f example000_simpleForward.png
```

Expected: `All checks passed!`; `N files already formatted`; ty: `All checks passed!`; `197 passed` (the `.pyx` files changed, hence the rebuild; results identical); the example runs to completion with exit status 0 (it writes its figure into the current directory, which the `rm` removes; the tracked figure in `examples/` is untouched). `git status --short` still shows only the 13 files and the config.

- [ ] **Step 6: Commit in two steps**

```bash
git add .pre-commit-config.yaml
git commit -m "chore: add pre-commit configuration"
git add examples/example000_simpleForward.py examples/example001_simpleBackward.py examples/example002_methodOrder.py examples/example003_noisyBackward.py examples/example004_fullComparison.py examples/example005_comparisonPyAbel.py examples/example006_simpleForwardAndBackward.py
git add src/openAbel/abel/base.pyx src/openAbel/abel/fmm.pyx src/openAbel/abel/hansenLaw.pyx src/openAbel/abel/trap.pyx src/openAbel/abel/wrap.pyx src/openAbel/mathFun.pyx
git commit -m "style: apply ruff and pre-commit formatting to the examples and Cython sources"
```

---

### Task 5: GitHub Actions CI and dependabot

**Files:**
- Create: `.github/workflows/ci.yml`, `.github/dependabot.yml`, `.github/workflows/dependabot-automerge.yml`
- Delete: `.travis.yml`

**Interfaces:**
- Consumes: `uv sync` (default groups `dev` + `docs`, Task 1), the lint commands of Task 4, `uv run mkdocs build --strict` (which only passes once Task 7 has landed; the workflow is exercised on GitHub, not here).
- Produces: a workflow named exactly `CI` (Task 6's `tag.yml.disabled` triggers on that name) at the path `.github/workflows/ci.yml` (the README badge links to it).

- [ ] **Step 1: Write `.github/workflows/ci.yml`**

```yaml
name: CI

on:
  pull_request:
    branches: [main]
  push:
    branches: [main]

jobs:
  lint:
    runs-on: ubuntu-latest
    env:
      UV_PYTHON: "3.14"
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v7
      - run: uv python install 3.14
      - run: uv sync
      - run: uv lock --check
      - run: uv run ruff check --no-fix
      - run: uv run ruff format --check
      - run: uv run ty check src tests
      - run: uv run mkdocs build --strict

  test:
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest]
        python-version: ["3.12", "3.13", "3.14", "3.14t"]
        include:
          - os: macos-latest
            python-version: "3.14"
    env:
      UV_PYTHON: ${{ matrix.python-version }}
    steps:
      - uses: actions/checkout@v7
      - uses: astral-sh/setup-uv@v7
      - run: uv python install ${{ matrix.python-version }}
      - run: uv sync
      - run: uv run pytest
```

(`--no-fix` matters: `pyproject.toml` sets `fix = true`, so a plain `ruff check` in CI would auto-fix and pass. `UV_PYTHON` pins the interpreter per leg so the free-threaded leg cannot pick up the runner's system Python.)

- [ ] **Step 2: Write `.github/dependabot.yml`**

```yaml
version: 2
updates:
  - package-ecosystem: "uv"
    directory: "/"
    schedule:
      interval: "weekly"
    groups:
      dev-dependencies:
        dependency-type: "development"
      production-dependencies:
        dependency-type: "production"
  - package-ecosystem: "github-actions"
    directory: "/"
    schedule:
      interval: "weekly"
    groups:
      actions:
        patterns: ["*"]
```

- [ ] **Step 3: Write `.github/workflows/dependabot-automerge.yml`**

```yaml
name: Dependabot auto-merge

on: pull_request

permissions:
  contents: write
  pull-requests: write

jobs:
  automerge:
    runs-on: ubuntu-latest
    if: github.actor == 'dependabot[bot]'
    steps:
      - run: gh pr merge --auto --squash "$PR_URL"
        env:
          PR_URL: ${{ github.event.pull_request.html_url }}
          GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}
```

- [ ] **Step 4: Delete the Travis configuration and validate the YAML**

```bash
git rm -q .travis.yml
uv run python -c "import pathlib, yaml; [yaml.safe_load(p.read_text()) for p in pathlib.Path('.github').rglob('*.yml')]; print('yaml ok')"
uv run pre-commit run --all-files
```

Expected: `yaml ok` (PyYAML comes with mkdocs); every hook `Passed` or `Skipped`; `git status --short` shows `D  .travis.yml` and the three untracked `.github` files, nothing else.

- [ ] **Step 5: Commit**

```bash
git add .github/workflows/ci.yml .github/dependabot.yml .github/workflows/dependabot-automerge.yml
git commit -m "ci: add GitHub Actions workflow and dependabot"
```

(`git rm` staged the deletion.)

---

### Task 6: Dormant publish, tag and docs workflows

**Files:**
- Create: `.github/workflows/publish.yml.disabled`, `.github/workflows/tag.yml.disabled`, `.github/workflows/docs.yml.disabled`

**Interfaces:**
- Consumes: the `CI` workflow name (Task 5), `[tool.cibuildwheel]` in `pyproject.toml` (Task 1), `mike` from the `docs` group (Task 1; `mkdocs.yml` lands in Task 7).
- Produces: nothing consumed later. The files are inert until renamed; this plan never renames them.

The `.disabled` suffix makes GitHub ignore the files. A commented-out `*.yml` is still parsed and shows up as a failing "invalid workflow" run, which is why disabling-by-rename is used. Dependabot does not see `.disabled` files either, hence the "refresh the pins" item in the checklist.

- [ ] **Step 1: Write `.github/workflows/publish.yml.disabled`**

```yaml
# DISABLED. This file ends in `.disabled` so GitHub ignores it. A commented-out `*.yml` is still parsed by
# GitHub and shows up as a failing "invalid workflow" run, which is why disabling-by-rename is used instead.
#
# Activation checklist (spec section 6):
#   1. Make the repository public (or accept a private release).
#   2. Register the Trusted Publisher on PyPI: project `openAbel`, owner `oliverhaas`, workflow
#      `publish.yml`, environment `pypi`.
#   3. Create the `pypi` environment in the GitHub repository settings.
#   4. Enable GitHub Pages from the `gh-pages` branch (docs.yml deploys there through mike).
#   5. Refresh the action pins in the three `.disabled` files: dependabot does not see them.
#   6. Rename `publish.yml.disabled`, `tag.yml.disabled` and `docs.yml.disabled` to `.yml`.
#   7. Delete the Read the Docs project or leave a redirect; set `[project.urls] Documentation` in pyproject.toml.
name: Publish

on:
  push:
    tags: ["v*"]
  workflow_dispatch:
    inputs:
      version:
        description: "Version to publish (without the v prefix)"
        required: true
        type: string

jobs:
  test:
    runs-on: ubuntu-latest
    env:
      UV_PYTHON: "3.14"
    steps:
      - uses: actions/checkout@v7
        with:
          ref: ${{ github.event.inputs.version && format('v{0}', github.event.inputs.version) || github.ref }}
      - uses: astral-sh/setup-uv@v7
      - run: uv python install 3.14
      - run: uv sync
      - run: uv run ruff check --no-fix
      - run: uv run ruff format --check
      - run: uv run ty check src tests
      - run: uv run pytest

  build_wheels:
    name: Wheels on ${{ matrix.os }}
    needs: test
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-latest]
    steps:
      - uses: actions/checkout@v7
        with:
          ref: ${{ github.event.inputs.version && format('v{0}', github.event.inputs.version) || github.ref }}
      - uses: astral-sh/setup-uv@v7
      # cibuildwheel >= 3.4 builds the cp314t targets of [tool.cibuildwheel] without an enable flag; v4 removed it.
      - uses: pypa/cibuildwheel@v3.4.1
      - uses: actions/upload-artifact@v7
        with:
          name: wheels-${{ matrix.os }}
          path: ./wheelhouse/*.whl

  build_sdist:
    needs: test
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v7
        with:
          ref: ${{ github.event.inputs.version && format('v{0}', github.event.inputs.version) || github.ref }}
      - uses: astral-sh/setup-uv@v7
      - run: uv build --sdist
      - uses: actions/upload-artifact@v7
        with:
          name: sdist
          path: ./dist/*.tar.gz

  publish:
    needs: [build_wheels, build_sdist]
    runs-on: ubuntu-latest
    environment:
      name: pypi
      url: https://pypi.org/project/openAbel/
    permissions:
      id-token: write
    steps:
      - uses: actions/download-artifact@v8
        with:
          path: dist
          pattern: "*"
          merge-multiple: true
      - uses: pypa/gh-action-pypi-publish@release/v1
```

- [ ] **Step 2: Write `.github/workflows/tag.yml.disabled`**

```yaml
# DISABLED. This file ends in `.disabled` so GitHub ignores it. A commented-out `*.yml` is still parsed by
# GitHub and shows up as a failing "invalid workflow" run, which is why disabling-by-rename is used instead.
# Activation checklist: see publish.yml.disabled.
#
# After a green CI run on main: if pyproject.toml's version is neither on PyPI nor tagged yet, tag it and
# dispatch the publish and docs workflows. Tag-exists is treated as "do nothing"; a broken tag must be deleted
# manually before a re-release. GITHUB_TOKEN-pushed tags don't trigger downstream workflows, hence the explicit
# dispatch.
name: Tag release

on:
  workflow_run:
    workflows: ["CI"]
    types: [completed]
    branches: [main]

permissions:
  contents: write
  actions: write

jobs:
  tag:
    if: github.event.workflow_run.conclusion == 'success' && github.event.workflow_run.event == 'push'
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v7
        with:
          fetch-depth: 2
      - uses: astral-sh/setup-uv@v7

      - name: Check version
        id: version
        run: |
          CURRENT_VERSION=$(uv version --short)
          echo "version=$CURRENT_VERSION" >> "$GITHUB_OUTPUT"
          if curl -sf "https://pypi.org/pypi/openAbel/$CURRENT_VERSION/json" > /dev/null; then
            echo "on_pypi=true" >> "$GITHUB_OUTPUT"
          else
            echo "on_pypi=false" >> "$GITHUB_OUTPUT"
          fi
          git fetch --tags
          if git rev-parse "v$CURRENT_VERSION" > /dev/null 2>&1; then
            echo "tag_exists=true" >> "$GITHUB_OUTPUT"
          else
            echo "tag_exists=false" >> "$GITHUB_OUTPUT"
          fi

      - name: Create and push tag
        if: steps.version.outputs.on_pypi == 'false' && steps.version.outputs.tag_exists == 'false'
        run: |
          git config user.name "github-actions[bot]"
          git config user.email "github-actions[bot]@users.noreply.github.com"
          git tag "v${{ steps.version.outputs.version }}"
          git push origin "v${{ steps.version.outputs.version }}"

      - name: Trigger publish and docs
        if: steps.version.outputs.on_pypi == 'false' && steps.version.outputs.tag_exists == 'false'
        env:
          GH_TOKEN: ${{ github.token }}
        run: |
          gh workflow run publish.yml -f version="${{ steps.version.outputs.version }}" -R ${{ github.repository }}
          gh workflow run docs.yml --ref "v${{ steps.version.outputs.version }}" -R ${{ github.repository }}
```

- [ ] **Step 3: Write `.github/workflows/docs.yml.disabled`**

```yaml
# DISABLED. This file ends in `.disabled` so GitHub ignores it. A commented-out `*.yml` is still parsed by
# GitHub and shows up as a failing "invalid workflow" run, which is why disabling-by-rename is used instead.
# Activation checklist: see publish.yml.disabled.
name: Docs

on:
  push:
    branches: [main]
    tags: ["v*"]
  workflow_dispatch:

permissions:
  contents: write
  pages: write
  id-token: write

jobs:
  deploy:
    runs-on: ubuntu-latest
    env:
      UV_PYTHON: "3.14"
    steps:
      - uses: actions/checkout@v7
        with:
          fetch-depth: 0
      - name: Configure Git
        run: |
          git config user.name "github-actions[bot]"
          git config user.email "github-actions[bot]@users.noreply.github.com"
      - uses: astral-sh/setup-uv@v7
      - run: uv python install 3.14
      - run: uv sync
      - name: Deploy docs (version tag)
        if: startsWith(github.ref, 'refs/tags/v')
        run: |
          VERSION=${GITHUB_REF#refs/tags/v}
          uv run mike deploy --push --update-aliases "$VERSION" latest
          uv run mike set-default --push latest
      - name: Deploy docs (main branch)
        if: github.ref == 'refs/heads/main'
        run: uv run mike deploy --push main
```

- [ ] **Step 4: Validate and commit**

```bash
uv run python -c "import pathlib, yaml; [yaml.safe_load(p.read_text()) for p in pathlib.Path('.github/workflows').glob('*.disabled')]; print('yaml ok')"
grep -L 'DISABLED' .github/workflows/*.disabled
uv run pre-commit run --all-files
git add .github/workflows/publish.yml.disabled .github/workflows/tag.yml.disabled .github/workflows/docs.yml.disabled
git commit -m "ci: add the dormant publish, tag and docs workflows"
```

Expected: `yaml ok`; `grep -L` prints nothing (every file carries the header); hooks clean; `git status --short` empty after the commit.

---

### Task 7: Documentation port to mkdocs-material

**Files:**
- Create: `mkdocs.yml`, `docs/index.md`, `docs/transform-types.md`, `docs/transform-methods.md`, `docs/remarks.md`, `docs/javascripts/mathjax.js`, `docs/examples/index.md`, `docs/examples/example000.md` ... `docs/examples/example005.md`, `docs/reference/api.md`, `docs/reference/changelog.md`
- Unchanged: the six `docs/examples/*.png`, `docs/superpowers/` (the Sphinx files were removed in Task 4)

**Interfaces:**
- Consumes: `README.md` (Task 1) and the formatted `examples/*.py` (Task 4), both included by `pymdownx.snippets` with `base_path: [".", "docs"]` and `check_paths: true`; `mkdocs`/`mkdocs-material`/`mike` from the `docs` group (Task 1).
- Produces: `uv run mkdocs build --strict` passing (the CI lint job of Task 5 runs it); the `site/` output directory (git-ignored).

Content rules: the pages below carry the text of the old RST pages with headings fixed and links pointed at the new page names; math uses `\(...\)` inline and `\[...\]` display (arithmatex generic mode, MathJax 3). Example pages include their script through a `{ .python }` fence with a `--8<--` snippet line. The changelog is the `0.7.0` entry for this branch.

- [ ] **Step 1: Confirm the starting state**

```bash
ls docs docs/examples
git status --short
```

Expected: `docs/` holds `examples/` and `superpowers/`; `docs/examples/` holds the six PNGs and nothing else; `git status --short` is empty. (If any `.rst` file or `docs/conf.py` is still present, Task 4 was not completed: stop and report.)

- [ ] **Step 2: Write `mkdocs.yml`**

```yaml
site_name: openAbel
site_description: Fast Abel transforms of equispaced data
site_url: https://oliverhaas.github.io/openAbel/
repo_url: https://github.com/oliverhaas/openAbel
repo_name: oliverhaas/openAbel
edit_uri: edit/main/docs/

# Internal planning notes, not part of the published site.
exclude_docs: |
  superpowers/

theme:
  name: material
  palette:
    - media: "(prefers-color-scheme: light)"
      scheme: default
      primary: red
      accent: red
      toggle:
        icon: material/brightness-7
        name: Switch to dark mode
    - media: "(prefers-color-scheme: dark)"
      scheme: slate
      primary: red
      accent: red
      toggle:
        icon: material/brightness-4
        name: Switch to light mode
  features:
    - navigation.instant
    - navigation.instant.progress
    - navigation.tracking
    - navigation.sections
    - navigation.expand
    - content.code.copy
    - content.code.annotate

markdown_extensions:
  - admonition
  - attr_list
  - md_in_html
  - pymdownx.arithmatex:
      generic: true
  - pymdownx.details
  - pymdownx.highlight:
      anchor_linenums: true
      line_spans: __span
      pygments_lang_class: true
  - pymdownx.inlinehilite
  - pymdownx.snippets:
      base_path: [".", "docs"]
      check_paths: true
  - pymdownx.superfences
  - toc:
      permalink: true

extra_javascript:
  - javascripts/mathjax.js
  - https://unpkg.com/mathjax@3/es5/tex-mml-chtml.js

extra:
  version:
    provider: mike

nav:
  - Home: index.md
  - Transform types: transform-types.md
  - Transform methods: transform-methods.md
  - Remarks: remarks.md
  - Examples:
      - examples/index.md
      - examples/example000.md
      - examples/example001.md
      - examples/example002.md
      - examples/example003.md
      - examples/example004.md
      - examples/example005.md
  - Reference:
      - API: reference/api.md
      - Changelog: reference/changelog.md
```

- [ ] **Step 3: Write `docs/javascripts/mathjax.js`**

```javascript
window.MathJax = {
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex",
  },
};

document$.subscribe(() => {
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
```

- [ ] **Step 4: Write `docs/index.md`**

The whole file is one line (the README is the single source of the landing page):

```markdown
--8<-- "README.md"
```

- [ ] **Step 5: Write `docs/transform-types.md`**

````markdown
# Transform types

In **openAbel** due to the equispaced discretization all methods truncate the Abel transform integral, e.g. for the
forward Abel transform

\[
F(y)=2\int_y^\infty\frac{f(r)r}{\sqrt{r^2-y^2}}dr\approx2\int_y^R\frac{f(r)r}{\sqrt{r^2-y^2}}dr\; .
\]

This is sometimes called finite Abel transform. Since \(f(r)\) often has compact support or decays very quickly (and
\(R\) can be chosen very large with a fast transform method) this introduces an arbitrarily small error.

Often one can use variable transformations or other discretizations to simplify the
calculation of the above integrals. However, often one is interested in exactly the in **openAbel** implemented case
on equispaced discretization. This is often due to the
[relation of the Abel transform with the Fourier and Hankel transforms](https://en.wikipedia.org/wiki/Abel_transform#Relationship_to_the_Fourier_and_Hankel_transforms)
and the desire to use the same discretization as the FFT or a discrete convolution, or just by the given data (e.g.
from experiments).

The type of transform can be chosen by setting the `forwardBackward` parameter:

```python
import openAbel

abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize)
```

The parameter `stepSize` is the grid spacing of the equidistant grid, `nData` the length of the data input array, and
`shift` is an offset of the samples to the symmetry axis and can usually be only 0 or 0.5 (input in units of
`stepSize`).

## Forward Abel transform

The forward Abel transform is defined as

\[
F(y)=2\int_y^\infty\frac{f(r)r}{\sqrt{r^2-y^2}}dr\approx2\int_y^R\frac{f(r)r}{\sqrt{r^2-y^2}}dr\; .
\]

The forward Abel transform is chosen by setting `forwardBackward=-1`.

## Backward (or inverse) Abel transform

The backward (or inverse) Abel transform is defined as

\[
f(r)=-\frac{1}{\pi}\int_r^\infty\frac{F'(y)}{\sqrt{y^2-r^2}}dy\approx-\frac{1}{\pi}\int_r^R\frac{F'(y)}{\sqrt{y^2-r^2}}dy\; .
\]

**openAbel** takes care of taking the derivative of the input data supplied by the user. The backward Abel transform is
chosen by setting `forwardBackward=1`.

## Backward (or inverse) Abel transform with derivative input

The backward (or inverse) Abel transform with derivative input is defined as

\[
f(r)=-\frac{1}{\pi}\int_r^\infty\frac{g(y)}{\sqrt{y^2-r^2}}dy\approx-\frac{1}{\pi}\int_r^R\frac{g(y)}{\sqrt{y^2-r^2}}dy\; .
\]

In contrast to the normal backward Abel transform, **openAbel** expects to get the derivative as input by the user.
The backward Abel transform with derivative input is chosen by setting `forwardBackward=2`.

## Modified forward Abel transform

What we call the modified forward Abel transform in **openAbel** is defined as the integral

\[
H(y)=2\int_y^\infty\frac{h(r)y^2}{r^2\sqrt{r^2-y^2}}dr\approx2\int_y^R\frac{h(r)y^2}{r^2\sqrt{r^2-y^2}}dr\; .
\]

I encountered this integral when a radial electric field of an atom (which has a \(1/r^2\) singularity we want to
integrate properly) is projected instead of a simpler function like with the normal forward Abel transform. One could
just use the parameter `shift = 0.5` instead to avoid the singularity of the electric field, but if one incorporates
the singularity in the actual integral the convergence is much better. I recommend writing similar methods if one
encounters other types of singularities in the Abel transform.

The modified forward Abel transform is chosen by setting `forwardBackward=-2`.
````

- [ ] **Step 6: Write `docs/transform-methods.md`**

````markdown
# Transform methods

In **openAbel** there are a couple of different algorithms for the calculation of the Abel transforms implemented,
although most of them are just for comparisons and it is recommended to only use the default method.

The main two obstacles when calculating the transforms numerically are the singularity at \(r=y\) and the dependence
of the result on \(y\), meaning computational complexity is quadratic \(O(N^2)\) if one naively integrates. The main
difference between the implemented transforms is how those two issues are treated.

When creating the Abel transform object the `method` keyword argument can be provided to choose different transform
methods:

```python
import openAbel

abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize, method=3, order=2)
```

The methods with end corrections can do the transformation in different orders of accuracy by setting the `order`
keyword argument; all other methods ignore `order`. Note when we talk about \(n\) order accuracy we usually mean
\((n+1/2)\) order accuracy due to the square root in the Abel transform kernel. For higher order methods the
transformed function has to be sufficiently smooth to achieve the full order of convergence, and in very extreme cases
the transform become unstable if high order is used on non-smooth functions. The length of the data vector `nData` we
denote as \(N\) in the math formulas.

Overall cases where a user should use anything other than `method = 3` (default) and `order = 2` (default) to
`order = 5` will be very rare. For a detailed comparison of the methods it is recommended to look at
[example004_fullComparison](examples/example004.md).

## Desingularized trapezoidal rule

```python
# order keyword argument is ignored (only first order implemented)
abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize, method=0)
```

The desingularized trapezoidal rule is probably the simplest practicable algorithm. It subtracts the singularity and
integrates it analytically, and numerically integrates the remaining desingularized term with the trapezoidal rule. In
the implementation this is done to first order, i.e. for the forward Abel transform this leads to

\[
F(y)=2\int_{y}^{R}\frac{(f(r)-f(y))r}{\sqrt{r^2-y^2}}dr+f(y)\sqrt{R^2-y^2}\;.
\]

Now the singularity seems to be removed, but a closer look and one can see that the singularity is still there in the
derivative of the integrand, so the convergence is first order in \(N\) instead of second order expected when using
trapezoidal rule. One can analytically remove the singularity in higher order with more terms, but this gets kinda
complicated (and possibly unstable, plus there are other practical issues). The trapezoidal rule portion of the method
leads to quadratic \(O(N^2)\) computational complexity of the method.

## Hansen-Law method

```python
# order keyword argument is ignored (only somewhat first order implemented)
abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize, method=1)
```

The Hansen-Law method by [Hansen and Law](https://www.osapublishing.org/josaa/abstract.cfm?uri=josaa-2-4-510) is a
space state model approximation of the Abel transform kernel. With that method recursively transforms a piecewise
linear approximation of the input functions to integrate analytically piece by piece. In principle this results in an
2nd order accurate transform, but the approximation of the Abel transform kernel as a sum of exponentials is quite
difficult. In other words the approximation

\[
\frac{1}{\sqrt{1-\exp{(-2t)}}}\approx\sum_{k=1}^K\exp{(-\lambda_kt)}
\]

is in practice not possible to achieve with high accuracy and reasonable \(K\). This is the main limitation of the
method, and the original space state model approximation has a typical relative error of \(10^{-3}\) at best -- then it
just stops converging with increasing \(N\). If one ignores several details that makes the method apparently linear
\(O(N)\) computational complexity, so it is implemented here for comparisons. More comments in the
[remarks](remarks.md).

## Trapezoidal rule with end corrections

```python
# 0 < order < 20
abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize, method=2, order=2)
```

The trapezoidal rule with end correction improves on the desingularized trapezoidal rule. It doesn't require
analytical integration because it uses precalculated end correction coefficients of arbitrary order. As described in
[Kapur](https://epubs.siam.org/doi/abs/10.1137/S0036142995287847) one can construct \(\alpha_i\) and \(\beta_i\) such
that the approximation

\[
\int_{a}^{b}f(x)dx \approx h\cdot\sum_{i=1}^{N-2}f(x_i) +
                           h\cdot\sum_{i=0}^{M-1}\alpha_if(x_{i-p}) +
                           h\cdot\sum_{i=0}^{M-1}\beta_if(x_{N-1-q})
\]

is accurate to order \(M\). Note that \(p\) and \(q\) should be chosen such that the correction is centered around the
end points: Similar to central finite differences this leads to an arbitrary order stable scheme, and thus incredibly
fast convergence and small errors. Otherwise it's not recommended to go higher than \(M=5\), again similar to forward
and backward finite differences. The trapezoidal rule portion of the method leads to quadratic \(O(N^2)\) computational
complexity of the method.

Since the calculation of the end correction coefficients requires some analytical calculations, is quite troublesome
and time consuming, they have been precalculated in *Mathematica* and stored in binary *\*.npy*, so they are only
loaded by the **openAbel** code when needed and don't have to be calculated. The
[*Mathematica* notebook](https://github.com/oliverhaas/openAbel/tree/main/add/calcEndCorr.nb) which was used to
calculate these end correction coefficients can be found in this repository as well.

## Fast Multipole Method with end corrections

```python
# 0 < order < 20
abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize, method=3, order=2)
```

The default and recommended method is the Fast Multipole Method (FMM) with end corrections. This method provides a
fast linear \(O(N)\) computational complexity transform of arbitrary order. The specific FMM used is based on Chebyshev
interpolation and nicely described and applied by
[Tausch](https://link.springer.com/chapter/10.1007/978-3-642-25670-7_6) on a similar problem. In principle the FMM
uses a hierarchic decomposition to combine a linear amount of direct short-range contributions and smooth
approximations of long-range contributions with efficient reuse of intermediate results to get in total a linear
\(O(N)\) computational complexity algorithm. This method thus provides extremely fast convergence and fast computation,
and is optimal in that sense for the intended purpose.

## Remarks on transforms of noisy data

For specifically the inverse Abel transform of noisy data (e.g. experimental data) there are a lot of algorithms
described in literature which might perform better in some aspects, since they either incorporate some assumptions
about the data or some kind of smoothing/filtering of the noise. A nice starting point for people interested in those
methods is the Python module [PyAbel](https://github.com/PyAbel/PyAbel).

However, there is no reason not to combine the methods provided in **openAbel** with some kind of filtering for nicer
results. I've had good results with [maximally flat filters](https://ieeexplore.ieee.org/document/7944698/), as seen
in [example003_noisyBackward](examples/example003.md), and with additional material in the
[Mathematica notebook](https://github.com/oliverhaas/openAbel/tree/main/add/calcMaxFlat.nb).

Overall even in this special case there are no algorithms to my knowledge which perform inherently better than the
default algorithms of **openAbel** by default.
````

- [ ] **Step 7: Write `docs/remarks.md`**

````markdown
# Remarks

In this section I want to mainly give some fairly informal remarks. Many are on different Abel transform algorithms
described in literature, and for those I'm giving the reasoning on why these algorithms were not useful to me. The
remaining remarks are just a collection of findings of me related to the Abel transform I thought worth mentioning.
Some of the remarks are fairly subjective for the use cases I was and am interested in. While this
section is meant to be fairly informal, I try to be direct, especially since I feel like many publications are
somewhat misleading, e.g. algorithms are often mislabeled as "fast" or "similar to the Fast Fourier Transform" even
though the computational complexity is nowhere near similar to the Fast Fourier Transform.

## PyAbel Python package

There is a Python package called [PyAbel](https://github.com/PyAbel/PyAbel), which focuses mainly on the inverse (or
backward/reconstruction/etc.) Abel transform. It's open source and pretty well documented and has several nice and
communicative developers which are active on the **PyAbel** GitHub repository.

If the reader is interested in testing many of the algorithms mentioned in this remarks section, I can only recommend
to look at **PyAbel**, as it does implement more algorithms for the inverse Abel transform than **openAbel**
(basically all except the main recommended one in **openAbel**). The
[example005_comparisonPyAbel](examples/example005.md) of **openAbel** uses **PyAbel** to compare many of the
algorithms talked about here as well.

## Hansen-Law method

Although the [Hansen-Law method](https://www.osapublishing.org/josaa/abstract.cfm?uri=josaa-2-4-510) is still
implemented in **openAbel**, it unfortunately has one major flaw (and one slightly problematic one). I think in
principle the space state model approximation is a pretty clever idea. That's why I actually tried to improve the
algorithm or apply the basic idea to other integrals (both to no success, yet). Unfortunately several things have to
fall into place for it to work nicely, which isn't the case for the Abel transform.

The first approximation which Hansen and Law make is -- roughly speaking -- that they implicitly approximate the data
by a piecewise polynomial, specifically piecewise linear in the original formulation. Going to higher order gets messy
quickly, but is in principle possible. I successfully tried that, but due to the next approximation (which is somewhat
flawed) it's not useful.

The second approximation is to rewrite the Abel transform kernel and approximate it by a sum of exponentials. At first
it looks like one could get a linear computational complexity \(O(N)\) algorithm. Problem is that the kernel has (even
after rewriting) a singularity, so it's obviously pretty difficult to approximate a singularity by a sum of
exponentials (again I tried many published approaches for that; most have flaws as well or are at least difficult and
don't lead to good enough results). Increasing the number of exponentials used increases the computational
complexity. My guess is that the algorithm is \(O(N \log(N))\) at best because of the increasing number of required
exponentials. In practice it turns out it's pretty much impossible to get anything really universally useful, unless
one aims only for fairly large errors (like Hansen-Law with roughly \(10^{-3}\)). For many use cases this is probably
enough (e.g. experimental ones which are usually fairly noisy anyway), so the method has still some value. But even if
that is the case and one ignores every problem I mentioned here the method is still not more efficient than the main
**openAbel** methods, so it's never the best choice as long as one does not have to implement the algorithms (FMM is a
lot more work to implement than Hansen-Law). And I have to mention that I did a lot of thinking on the topic, and for
example one could subtract the singularity in the impulse response, but this usually leads to not useful algorithms
or the use of basically the same approaches and algorithms as **openAbel** to make it competitive, just one is taking
a huge detour.

## Desingularized quadrature

Like the Hansen-Law method the desingularized quadrature is implemented in **openAbel**, but still mostly a remnant
of testing different Abel transform methods and I don't recommend using it. The basic idea of desingularizing
integrals should be familiar to almost anyone who tried to numerically integrate a function with a singularity and
wanted to improve convergence. In context of a problem similar to the Abel transform it is discussed by
[Tausch](https://link.springer.com/chapter/10.1007/978-3-642-25670-7_6), which is one of the main references for the
Fast Multipole Method in **openAbel** as well. I actually tried several higher order versions of this, and the effort
is not really worth it, since the end corrections used in **openAbel** are much more efficient. And it gets very
complicated -- maybe impossible -- to program if one tries to avoid instabilities; I'm actually not sure if my test
implementations were definitely reliable. And of course this method, it's \(O(N^2)\), is slower than the main
**openAbel** methods.

## Piecewise polynomial analytic integration

There are many different publications which basically use the same idea: Interpolate the data by a piecewise
polynomial and use the known analytic integral for every polynomial piece.
[Dasch](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-31-8-1146) calls it onion peeling or filtered back
projection, [Bordas](https://aip.scitation.org/doi/abs/10.1063/1.1147044) does basically the same, and there are
probably more publications. I actually thought initially when I decided to use the Fast Multipole Method, that the
piecewise polynomial analytic integration would be useful in combination. In principle it works, but is again pretty
messy and the end corrections used in **openAbel** are overall much more efficient in every sense once implemented.
But without the Fast Multipole Method it's a slow \(O(N^2)\) method as well, and that slow way is what all
publications do to my knowledge.

## Analytic integration of a basis set expansion

There are many publications which use some kind of basis set expansion applied to the data, then use analytic
transform of each basis function to construct the total transform. So this is similar to **Piecewise polynomial
analytic integration**, but with a basis for the whole domain and not just piecewise.

Often a polynomial basis set expansion for the whole data set and then transform each basis polynomial analytically.
This obviously only works well (regarding error) if the data has somewhat polynomial behavior regarding the whole
domain. Since one can choose orthogonal polynomial basis sets the expansion is at least fairly fast, but since the
analytical transforms of polynomials are not very nice this approach overall is not very efficient. In a specific case
where the data is basically a low order polynomial this of course would work really well, but in the general case it's
not useful.

Other basis sets might have nicer transforms, but are not orthogonal, so the expansion of the data is more difficult.
I tried to find some "good" basis, but in one way or another one shifts the difficulty to another area, e.g. function
approximation, and I did not get a useful approach. Again, for some very specific data sets one might find a very
small but usable basis set.

One example of such an algorithm in literature is the
[BASEX algorithm by Dribinski](https://aip.scitation.org/doi/abs/10.1063/1.1482156). In this method a "Gaussian" basis
set -- just to note it's somewhat Gaussian, not the "normal" Gaussian -- is used. It seems to be very popular, as the
publication has 791 citations as of writing this. I'm guessing mainly because the code was freely available and the
basis set implicitly applied some smoothing in the transform, which usually produces nicer pictures without tweaking
than other algorithms. I'm fairly convinced that one can achieve similarly nice results with other methods and some
smoothing. Similar to other methods described here the method can be tested in **PyAbel**. The preprocessing is
incredibly painfully slow (it's \(O(N^3)\) I think, and it takes minutes for even small arrays \(N=1000\), where
**openAbel**'s main methods are \(O(N)\) and take milliseconds), and the actual transform is not much better
(\(O(N^2)\) and **openAbel** is \(O(N)\) again). Overall **BASEX** is a fairly often cited algorithm nevertheless.

I can see how in some cases one might be able to choose a nicely suitable basis set to enforce some structure in either
the projected or reconstructed data. I expect this would be the only case where such an approach would make sense,
but this is very problem specific and thus much less universal than the main methods of **openAbel** intend to be.
One example of such an approach is described by [Gerber](https://aip.scitation.org/doi/10.1063/1.4793404), and often
called linBASEX (e.g. in **PyAbel**). Due to the underlying physical process Gerber expects or knows that his data has
some structure, and enforces it by choosing a specific basis set. In the general case one could probably achieve
similar results by fitting the expected structure basis set to the data and then using accurate black-box Abel
transform functions like in **openAbel**.
````

- [ ] **Step 8: Write `docs/examples/index.md`**

```markdown
# Examples

The scripts below live in the
[`examples/`](https://github.com/oliverhaas/openAbel/tree/main/examples) directory of the repository and are
reproduced here with their output figures. They need `matplotlib`; example005 also needs `PyAbel`.

- [example000_simpleForward](example000.md): a forward transform of a Gaussian.
- [example001_simpleBackward](example001.md): a backward transform of a Gaussian, with and without analytic
  derivative input.
- [example002_methodOrder](example002.md): switching transform methods and orders.
- [example003_noisyBackward](example003.md): filtering and transforming noisy data.
- [example004_fullComparison](example004.md): accuracy and timing comparison of all **openAbel** methods.
- [example005_comparisonPyAbel](example005.md): comparison of **openAbel** with **PyAbel** methods.
```

- [ ] **Step 9: Write the six example pages**

`docs/examples/example000.md`:

````markdown
# example000_simpleForward

This example is just a simple forward transform of a Gaussian. Aside from showing how to do a simple forward
transform, this example shows how for non truncated domain an error is introduced.

![Simple forward transform of a Gaussian.](example000_simpleForward.png)

```{ .python }
--8<-- "examples/example000_simpleForward.py"
```
````

`docs/examples/example001.md`:

````markdown
# example001_simpleBackward

This example is just a simple backward transform of a Gaussian. Aside from showing how to do a simple backward
transform, this example shows how for non truncated domain an error is introduced, and how taking the derivative
analytically of the data (if possible) improves the error. Taking numerical derivatives always amplifies noise and
increases the resulting error.

![Simple backward transform of a Gaussian.](example001_simpleBackward.png)

```{ .python }
--8<-- "examples/example001_simpleBackward.py"
```
````

`docs/examples/example002.md`:

````markdown
# example002_methodOrder

This example shows how to switch to other transform methods and orders. It illustrates how quickly the errors of high
order methods converge to machine precision, even for very small data sets. It is of course important that the input
data is sufficiently smooth and other errors (e.g. truncation errors) are small enough as well.

![Different methods and orders.](example002_methodOrder.png)

```{ .python }
--8<-- "examples/example002_methodOrder.py"
```
````

`docs/examples/example003.md`:

````markdown
# example003_noisyBackward

This example shows how to filter and transform noisy data. It illustrates how noisy input data can lead to large
errors of the backward transform result, and how filters can be used to -- at least visually -- alleviate those
errors.

![Backward transform of noisy data.](example003_noisyBackward.png)

Maximally flat filters have been calculated as described in a paper by
[Hosseini](https://ieeexplore.ieee.org/document/7944698/), and a small
[Mathematica script](https://github.com/oliverhaas/openAbel/blob/main/add/calcMaxFlat.nb) of the calculation is
provided in the additional materials.

```{ .python }
--8<-- "examples/example003_noisyBackward.py"
```
````

`docs/examples/example004.md`:

````markdown
# example004_fullComparison

This example provides a rather extensive comparison of different **openAbel** methods. It shows how well the main
methods perform in every regard, especially if the input data is sufficiently smooth. Most importantly one can see
the fast convergence of the higher order methods in the bottom left plot, and the linear computational complexity of
the main methods in the bottom middle and right plots.

![Comparison of different openAbel methods.](example004_fullComparison.png)

```{ .python }
--8<-- "examples/example004_fullComparison.py"
```
````

`docs/examples/example005.md`:

````markdown
# example005_comparisonPyAbel

This example provides a rather extensive comparison of different **openAbel** methods with
[PyAbel](https://github.com/PyAbel/PyAbel) methods. It shows how well the main methods of **openAbel** perform in
every regard, especially if the input data is sufficiently smooth.

Since **PyAbel**'s focus is on the backward (or inverse) transform, this example does it as well.

![Comparison of different openAbel with PyAbel methods.](example005_comparisonPyAbel.png)

```{ .python }
--8<-- "examples/example005_comparisonPyAbel.py"
```
````

- [ ] **Step 10: Write `docs/reference/api.md`**

````markdown
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

Raises `ValueError` if a parameter has a non-viable value (for example `order <= 0`, or too few data points for the
requested order) and `NotImplementedError` if the chosen method does not support the given parameters (for example
an unknown `method`, a `shift` other than `0.0` or `0.5` with methods `2` and `3`, or the modified forward transform
with the Hansen-Law method).

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
Samples beyond that are ignored. For the Hansen-Law method (`method=1`) the boundary arguments are ignored.

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
````

- [ ] **Step 11: Write `docs/reference/changelog.md`**

```markdown
# Changelog

## 0.7.0 (unreleased)

### Packaging and tooling

- Python >= 3.12; wheels for CPython 3.12, 3.13, 3.14 and the free-threaded 3.14t on Linux x86_64 and macOS arm64.
- `pyproject.toml` (PEP 621) with a `src/` layout; `setup.py` only compiles the Cython extensions. Version `0.7.0` is
  exposed as `openAbel.__version__`.
- Cython 3, numpy 2 and scipy >= 1.13.
- pytest replaces nose; ruff, ty and pre-commit; GitHub Actions CI (Linux 3.12-3.14t, macOS 3.14) replaces Travis.
- Documentation ported from Sphinx to mkdocs-material.

### Fixed

- Every `Abel(...)` construction died with exit status 255 on glibc >= 2.38 (e.g. Ubuntu 24.04): the internal
  allocator requested alignment 0 from `aligned_alloc` and then called `exit` on the NULL it got back. The allocator
  now uses 64-byte alignment, rounds the size up as C11 requires, and raises `MemoryError` on failure.
- Backward transform with `method=0` raised `FileNotFoundError` (wrong coefficient path).
- Modified forward transform (`forwardBackward=-2`) with `shift=0.5` raised `KeyError` with `method=2` and crashed
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

### Changed

- Cython 3 build: exception clauses moved after `nogil`, `cpow=True` keeps the integer power semantics of the FMM
  code, `freethreading_compatible=True`.
- The backward transform with `method=0` was never usable before this release; its first-order accuracy is
  low (relative error around 7e-2 on the Gaussian test case), which is expected for the method.

## Earlier versions

No changelog was kept before 0.7.0.
```

- [ ] **Step 12: Build the site strictly and check the includes**

```bash
uv run mkdocs build --strict
grep -c 'import openAbel' site/examples/example000/index.html
grep -c 'arithmatex' site/transform-types/index.html
grep -c 'Quick start' site/index.html
ls site | grep -c superpowers
uv run pre-commit run --all-files
```

Expected: the build ends with `Documentation built in ... seconds` and exit status 0. (mkdocs-material prints a boxed "Warning from the Material for MkDocs team" about MkDocs 2.0 first; that banner is informational, `--strict` only fails on MkDocs warnings, and a failure shows as `Aborted with N warnings in strict mode`.) Then a count of at least `1` three times (the example script, the math markup and the README landed in the pages), `0` (the planning notes are excluded from the site), and every hook clean. `site/` is git-ignored.

- [ ] **Step 13: Commit**

```bash
git add mkdocs.yml docs/index.md docs/transform-types.md docs/transform-methods.md docs/remarks.md docs/javascripts/mathjax.js
git add docs/examples/index.md docs/examples/example000.md docs/examples/example001.md docs/examples/example002.md docs/examples/example003.md docs/examples/example004.md docs/examples/example005.md
git add docs/reference/api.md docs/reference/changelog.md
git commit -m "docs: port the documentation to mkdocs-material"
```

`git status --short` is empty afterwards.

---

### Task 8: Final verification

No new files. Every check below must pass on the branch as it stands; anything that needs a fix is a defect of an earlier task, to be fixed in a follow-up commit here (Conventional Commits, named files) and reported.

**Files:**
- Read-only over the whole tree; `dist/` is produced (git-ignored).

**Interfaces:**
- Consumes: everything from Tasks 1-7.
- Produces: the verified branch `feat/modern-tooling`, ready for the merge decision.

- [ ] **Step 1: Hooks, lint, types, docs**

```bash
uv run pre-commit run --all-files
uv run ruff check --no-fix
uv run ruff format --check
uv run ty check src tests
uv run mkdocs build --strict
uv lock --check
git status --short
```

Expected: every hook `Passed`/`Skipped`; `All checks passed!` twice; formatted; `Documentation built in ...`; the lock is up to date; `git status --short` prints nothing.

- [ ] **Step 2: The suite on every interpreter**

```bash
uv sync --python 3.12 --reinstall-package openAbel && uv run pytest -q
uv sync --python 3.13 --reinstall-package openAbel && uv run pytest -q
uv sync --python 3.14t --reinstall-package openAbel && uv run pytest -q
uv sync --python 3.14 --reinstall-package openAbel && uv run pytest -q
```

Expected: `197 passed` four times; `.venv` ends on 3.14.

- [ ] **Step 3: Build the distributions**

```bash
rm -rf dist
uv build
ls dist
tar tzf dist/openabel-0.7.0.tar.gz | grep -c '\.npy$'
tar tzf dist/openabel-0.7.0.tar.gz | grep -c '\.pyx$'
unzip -l dist/openabel-0.7.0-*.whl | grep -c '\.npy$'
unzip -l dist/openabel-0.7.0-*.whl | grep -c '\.so$'
unzip -l dist/openabel-0.7.0-*.whl | grep -c '\.c$'
```

Expected: `dist/openabel-0.7.0.tar.gz` and one wheel `dist/openabel-0.7.0-cp314-cp314-linux_x86_64.whl` (the tag matches the interpreter uv built with; `cp313` is also acceptable if uv picked 3.13). Counts: `238` `.npy` in the sdist, `8` `.pyx` in the sdist, `238` `.npy` in the wheel, `8` `.so` in the wheel, `0` `.c` in the wheel.

- [ ] **Step 4: Install the sdist into a fresh environment and run the suite against it**

```bash
uv run --no-project --isolated --python 3.14 --with dist/openabel-0.7.0.tar.gz --with pytest python -m pytest tests -q -p no:cacheprovider
```

Expected: the sdist is compiled in an isolated environment (a few minutes) and the run ends with `197 passed in ...`. This is the check that the 2020 sdist failed (it shipped `.c` files that `setup.py` did not reference).

- [ ] **Step 5: Branch summary**

```bash
git log --oneline main..HEAD
git status --short
```

Expected: the commits of Tasks 1-7 in order (plus the spec and plan commits below them), nothing uncommitted. Report the log and every measured value (test counts, artifact names, `.npy` counts) in the task report.
