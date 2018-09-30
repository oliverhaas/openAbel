
openAbel README
=========

**Note:** It's best to view this readme in the 
`openAbel Documentation <https://openabel.readthedocs.io/en/latest/index.html>`_




Introduction
--------------


The main goal of **openAbel** is to provide fast and efficients Abel transforms of equispaced data
from smooth functions in Python with all actual calculations done in Cython. 
The most useful methods implemented in this module for that purpose use the Fast Multipole Method combined with
arbitrary order end correction of the trapezoidal rule to achieve small errors and fast convergence,
as well as linear computational complexity. A couple of other methods are implemented for comparisons.
Abel transform function can be used from Python and with numpy arrays or from Cython using pointers.


..\u0020Table of contents
..\u0020--------------

..\u0020- `Quick Start`_
..\u0020- `Dependencies`_
..\u0020- `Issues`_
..\u0020- `Transform Methods`_
..\u0020- `Copyright and License`_


Quick Start
--------------

In most cases this should be pretty simple:

- Clone the repository: :code:`git clone https://github.com/oliverhaas/openAbel.git`
- Install: :code:`sudo python setup.py install`
- Run example: :code:`python examples/example001_Gaussian.py`

This assumes dependencies are already met and you run a more or less similar system to mine (see `Dependencies`_).


Dependencies
--------------

The code was run on several Ubuntu systems without problems. More specific I'm running Ubuntu 16.04 and the following libraries and
Python modules, which were all installed the standard way with either :code:`sudo apt install moduleName` or `sudo pip install moduleName`. 

- Python 2.7.12

- Numpy 1.14.0

- Scipy 1.1.0

- Cython 0.23.4

- Matplotlib 2.2.0

- gcc 5.4.0


As usual newer versions of the libraries should work as well, and many older versions will too. I'm sure it's possible to
get **openAbel** to run on vastly different systems, like e.g. Windows systems, but obviously I haven't extensively tested
different setups.


Issues
--------------

If there are any issues, bugs or feature request just let me know. As of now there are some gaps in the implementation, e.g.
not all transform types are available in all methods, but since the default method vastly outperforms every other method 
anyway it's not really a pressing issue.


..\u0020Transform Methods
..\u0020--------------

..\u0020The forward and backward (or inverse) Abel transforms are usually defined as

..\u0020.. image:: https://latex.codecogs.com/gif.latex?F(y)=2\int_y^\infty\frac{f(r)r}{\sqrt{r^2-y^2}}dr
..\u0020.. image:: https://latex.codecogs.com/gif.latex?f(r)=-\frac{1}{\pi}\int_r^\infty\frac{F'(y)}{\sqrt{y^2-r^2}}dy\;.
..\u0020    
..\u0020In this section we will mostly talk about the forward Abel transform and then give some remarks on the inverse Abel transform at the end.

..\u0020It should be noted that often one can use variable transformations or other discretizations to simplify the calculation of the above
..\u0020integrals. However, often one is interested in exactly the in **openAbel** implemented case on equispaced discretization. This is often due to the 
..\u0020`relation of the Abel transform with the Fourier and Hankel transforms <https://en.wikipedia.org/wiki/Abel_transform#Relationship_to_the_Fourier_and_Hankel_transforms>`_ and the desire to use the same discretization as the FFT or a discrete convolution, or just by the given data (e.g. from experiments).

..\u0020When solved numerically on equispaced data the input data *f(r)* has to be truncated, such that the upper bound of the integral is some
..\u0020finite value *R>0*. Usually this isn't a problem since in most cases the functions decay rapidly with *r*. We obtain the sometimes called finite Abel transform

..\u0020.. image:: https://latex.codecogs.com/gif.latex?F(y)=2\int_y^R\frac{f(r)r}{\sqrt{r^2-y^2}}dr\; .

..\u0020There are two obstacles when calculating the transforms numerically: If one wants the output data on the same grid as the input data on *N*
..\u0020grid points, computational complexity is quadratic *O(N^2)*. And the singularity at *r=y* is difficult to handle efficiently.

..\u0020The first problem can be solved by the `Fast Multipole Method <https://en.wikipedia.org/wiki/Fast_multipole_method>`_ (FMM). The main reference
..\u0020for the implementation done here was the description of the Chebyshev Interpolation FMM by `Tausch <https://link.springer.com/chapter/10.1007/978-3-642-25670-7_6>`_.
..\u0020This leads to an *O(N)* algorithm when applied to the (discretized and truncated) Abel transform.

..\u0020The second problem can be solved by end correction to the trapezoidal rule. This is somewhat related to the fairly widely known `Euler-Maclaurin formula <https://en.wikipedia.org/wiki/Euler%E2%80%93Maclaurin_formula>`_ or `Gregory rules <https://www.sciencedirect.com/science/article/pii/0377042794902968>`_. 
..\u0020End corrections for singular function are described in several publications, but the main reference of **openAbel** was a paper by `Kapur <https://epubs.siam.org/doi/abs/10.1137/S0036142995287847>`_.
..\u0020If data points outside of the integration interval can be provided these end corrections are arbitrary order stable. Otherwise I wouldn't
..\u0020recommend going higher than 5th order. As of now we provide the coefficients up to 20th order.
..\u0020Since the calculation of the end correction coefficients requires some analytical calculations, is quite troublesome and time consuming, 
..\u0020they have been precalculated in *Mathematica* and stored in binary *\*.npy* , so they are only loaded by the **openAbel** code
..\u0020when needed and don't have to be calculated. The *Mathematica* `notebook <https://github.com/oliverhaas/openAbel/tree/master/add/calcEndCorr.nb>`_ which was used to calculate these end correction coefficients can be found in this repository as well.

..\u0020Overall to my knowledge there are no better methods for the described purpose.
..\u0020For specifically the inverse Abel transform of noisy data there are a lot of algorithms described in literature which perform better in
..\u0020some aspects, since they either incorporate some assumptions about the data or some kind of smoothing/filtering of the noise. A nice
..\u0020starting point for people interested in that is the Python module [PyAbel](https://github.com/PyAbel/PyAbel). However, one can use 
..\u0020**openAbel** for a noisy inverse transform as well, but one should do some manual filtering beforehand. I've had good results with
..\u0020`maximally flat filters <https://ieeexplore.ieee.org/document/7944698/>_ (see `notebook <https://github.com/oliverhaas/openAbel/tree/master/add/calcEndCorr.nb>`_[example003](examples/example003_inverse.py) 
..\u0020and the *Mathematica* [notebook](add/calcMaxFlat.nb)).

..\u0020In the examples some more details are discussed and mentioned; in general the examples are a good way to learn how to understand and
..\u0020use the code.


Copyright and License
--------------

Copyright &copy; 2016-2018 Oliver Sebastian Haas.

The code **openAbel** is published under the GNU GPL version 3. This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

For more information see the GNU General Public License copy provided in this repository `LICENSE <https://github.com/oliverhaas/openAbel/tree/master/LICENSE>`_.












