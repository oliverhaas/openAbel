# Examples

The scripts below live in the
[`examples/`](https://github.com/oliverhaas/openAbel/tree/main/examples) directory of the repository and are
reproduced here with their output figures. They need `matplotlib`; example005 also needs `PyAbel`. In a checkout the
`examples` dependency group provides both: `uv run --group examples python examples/example000_simple_forward.py`.

- [example000_simple_forward](example000.md): a forward transform of a Gaussian.
- [example001_simple_backward](example001.md): a backward transform of a Gaussian, with and without analytic
  derivative input.
- [example002_method_order](example002.md): switching transform methods and orders.
- [example003_noisy_backward](example003.md): filtering and transforming noisy data.
- [example004_full_comparison](example004.md): accuracy and timing comparison of all **openAbel** methods.
- [example005_comparison_pyabel](example005.md): comparison of **openAbel** with **PyAbel** methods.
