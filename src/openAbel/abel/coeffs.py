"""Precomputed end-correction and filter coefficients, loaded eagerly from ``coeffsData/*.npy`` at import."""

from pathlib import Path
from types import MappingProxyType

import numpy as np

dataDir = Path(__file__).parent / "coeffsData"


def loadCoeffs() -> MappingProxyType[str, MappingProxyType[int, np.ndarray]]:
    """Load every ``coeffsData/*.npy`` file into a read-only family -> order -> coefficients mapping."""
    # Outer key: coefficient family, i.e. the file name without its "_NN" suffix. Inner key: the order NN.
    loaded: dict[str, dict[int, np.ndarray]] = {}
    for path in sorted(dataDir.glob("*.npy")):
        coeffsName, _, order = path.stem.rpartition("_")
        loaded.setdefault(coeffsName, {})[int(order)] = np.load(path).astype(np.double)
    return MappingProxyType({name: MappingProxyType(byOrder) for name, byOrder in loaded.items()})


coeffsAllDict = loadCoeffs()


def getCoeffs(coeffsName: str, order: int) -> np.ndarray:
    """Return the coefficients of family ``coeffsName`` for the given ``order``."""
    return coeffsAllDict[coeffsName][order]
