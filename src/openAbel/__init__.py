"""openAbel: fast Abel transforms of equispaced data."""

from importlib.metadata import version

from .abel import Abel

__all__ = ["Abel"]
__version__ = version("openAbel")
