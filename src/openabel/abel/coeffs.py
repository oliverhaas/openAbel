"""Precomputed end-correction and filter coefficients, loaded eagerly from ``coeffs_data/*.npy`` at import."""

from pathlib import Path
from types import MappingProxyType

import numpy as np

DATA_DIR = Path(__file__).parent / "coeffs_data"


def load_coeffs() -> MappingProxyType[str, MappingProxyType[int, np.ndarray]]:
    """Load every ``coeffs_data/*.npy`` file into a read-only family -> order -> coefficients mapping."""
    # Outer key: coefficient family, i.e. the file name without its "_NN" suffix. Inner key: the order NN.
    loaded: dict[str, dict[int, np.ndarray]] = {}
    for path in sorted(DATA_DIR.glob("*.npy")):
        coeffs_name, _, order = path.stem.rpartition("_")
        loaded.setdefault(coeffs_name, {})[int(order)] = np.load(path).astype(np.double)
    return MappingProxyType({name: MappingProxyType(by_order) for name, by_order in loaded.items()})


COEFFS_ALL = load_coeffs()


def get_coeffs(coeffs_name: str, order: int) -> np.ndarray:
    """Return the coefficients of family ``coeffs_name`` for the given ``order``."""
    by_order = COEFFS_ALL.get(coeffs_name)
    if by_order is None or order not in by_order:
        raise ValueError(f"No {coeffs_name} coefficients of order {order}.")
    return by_order[order]
