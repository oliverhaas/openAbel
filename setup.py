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
