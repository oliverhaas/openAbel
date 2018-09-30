Transform Methods
=================


In **openAbel** there are a couple of different algorithms for the calculation of the Abel transforms implemented, although most of them are just for comparisons and is it recommended to only use the default method. 
Due to the equispaced discretization all methods truncate the Abel transform integral, e.g. for the forward Abel transform 

.. math::
        F(y)=2\int_y^\infty\frac{f(r)r}{\sqrt{r^2-y^2}}dr\approx2\int_y^R\frac{f(r)r}{\sqrt{r^2-y^2}}dr\; .
        
This is sometimes called finite Abel transform. Since :math:`f(r)` often has compact support or decays very quickly this usually introduces only a small error.

When creating the Abel transform object the :code:`method` keyword argument can be provided to chose different transform methods:

.. code-block:: python

    import openAbel
    abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize, method = 3, order = 2)

The methods with end corrections can do the transformation in different orders of accuracy by setting
the :code:`order` keyword argument; all other methods ignore :code:`order`. Note when we
talk about :math:`n` th order accuracy we usually mean :math:`(n+1/2)` th order accuracy due to the square
root in the Abel transform kernel. For higher order methods the transformed function has to be sufficiently
smooth to achieve the full order of convergence, and in very extreme cases the transform become unstable if
high order is used on non-smooth functions.
The length of the data vector :code:`nData` we denote as :math:`N` in the math formulas.


The different methods are (numbering coincides with the :code:`method` keyword argument integer):    

    0. The **desingularized trapezoidal rule** is probably the simplest practicable algorithm. 
    It subtracts the singularity and integrates it analytically, and numerically integrates the 
    remaining desingularized term with the trapezoidal rule. In the implementation this is done to first order, i.e.
    
    .. math::
            F(y)=2\int_{y}^{R}\frac{(f(r)-f(y))r}{\sqrt{r^2-y^2}}dr+f(y)\sqrt{R^2-y^2}\;.
            
    Now the singularity seems to be removed, but a closer look and one can see that the singularity
    is still there in the derivative of the integrand, so the convergence is first order in :math:`N`
    instead of second order expected when using trapezoidal rule. One can analytically remove the
    singularity in higher order with more terms, but for higher order this gets kinda complicated 
    (and possibly unstable, plus there are other issues). The trapezoidal rule portion of the method 
    leads to quadratic :math:`O(N^2)` computational complexity of the method.

    1. The **Hansen-Law method** by `Hansen and Law <https://www.osapublishing.org/josaa/abstract.cfm?uri=josaa-2-4-510>`_ 
    is a space state model approximation of the Abel transform kernel.
    With that method recursively transforms a piecewise linear approximation of the input functions 
    to integrate analytically piece by piece. In principle this results in an 2nd order accurate
    transform, but the approximation of the Abel transform kernel as a sum of exponentials is quite difficult.
    In other words the approximation

    .. math::
            \frac{1}{\sqrt{1-\exp{(-2t)}}}\approx\sum_{k=1}^K\exp{(-\lambda_kt)} 
        
    is in practice not possible to achieve with high accuracy and reasonable :math:`K`. This is the main 
    limitation of the method, and the original space state model approxmation has a typical relative
    error of :math:`1.e-3` at best -- then it just stops converging with increasing :math:`N`.
        
    2. The **trapezoidal rule with end corrections** improves on the desingularized trapezoidal rule.
    It doesn't require analytical integration because it uses precalculated end correction coefficients
    of arbitrary order. As described in `Kapur <https://epubs.siam.org/doi/abs/10.1137/S0036142995287847>`_
    one can contruct :math:`\alpha_i` and :math:`\beta_i` such that the approxmation

    .. math::
            \int_{a}^{b}f(x)dx \approx h\cdot\sum_{i=1}^{N-2}f(x_i) + 
                                       h\cdot\sum_{i=0}^{M-1}\alpha_if(x_{i-p}) + 
                                       h\cdot\sum_{i=0}^{M-1}\beta_if(x_{N-1-q})
    
    is accurate to order :math:`M`. Note that :math:`p` and :math:`q` should be chosen such that the correction is
    centered around the end points: Similar to central finite differences this leads to an arbitrary order stable scheme.
    Otherwise it's not recommended to go higher than :math:`M=5`, again similar to forward and backward finite
    differences. The trapezoidal rule portion of the method leads to quadratic :math:`O(N^2)` computational
    complexity of the method.
    
    3. The *default method* is the **Fast Multipole Method (FMM) with end corrections**. This method provides a fast
    linear :math:`O(N)` computational complexity transform of arbitrary order.
    The specific FMM used is based on Chebyshev interpolation and nicely described
    and applied by `Tausch <https://link.springer.com/chapter/10.1007/978-3-642-25670-7_6>`_ on a similar problem.
    In principle the FMM uses a hierarchic decomposition to combine a linear amount of direct short-range contributions
    and smooth approximations of long-range contributions with efficient reuse of intermediate results to get in total 
    a linear :math:`O(N)` computational complexity algorithm.



Overall cases where a user should use anything other than :code:`method = 3` and :code:`order = 2` to :code:`order = 5`
will be very rare.

For a detailed comparison of the methods it is recommended to look at exampleXXX.

