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
