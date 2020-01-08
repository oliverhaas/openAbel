Transform Types
=================


In **openAbel** due to the equispaced discretization all methods truncate the Abel transform integral, e.g. for the forward Abel transform 

.. math::
        F(y)=2\int_y^\infty\frac{f(r)r}{\sqrt{r^2-y^2}}dr\approx2\int_y^R\frac{f(r)r}{\sqrt{r^2-y^2}}dr\; .
        

This is sometimes called finite Abel transform. Since :math:`f(r)` often has compact support or decays very quickly (and :math:`R` can
be chosen very large with a fast transform method) this introduces an arbitrarily small error.

It should be noted that often one can use variable transformations or other discretizations to simplify the calculation of the above
integrals. However, often one is interested in exactly the in **openAbel** implemented case on equispaced discretization. This is often due to the 
`relation of the Abel transform with the Fourier and Hankel transforms <https://en.wikipedia.org/wiki/Abel_transform#Relationship_to_the_Fourier_and_Hankel_transforms>`_ and the desire to use the same discretization as the FFT or a discrete convolution, or just by the given data (e.g. from experiments).

The type of transform can be chosen by setting the :code:`forwardBackward` parameter:

.. code-block:: python

    import openAbel
    abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize)

The parameter :code:`stepSize` is the grid spacing of the equidistant grid, :code:`nData` the length of the data input array, and :code:`shift` is an offset of the samples to the symmetry axis and can usually be only 0 or 0.5 (input in units of :code:`stepSize`).


Forward Abel Transform
--------------
The Forward Abel transform is defined as

.. math::
        F(y)=2\int_y^\infty\frac{f(r)r}{\sqrt{r^2-y^2}}dr\approx2\int_y^R\frac{f(r)r}{\sqrt{r^2-y^2}}dr\; .

The Forward Abel Transform is chosen by setting :code:`forwardBackward=-1`.



Backward (or Inverse) Abel Transform
--------------
The Backward (or Inverse) Abel transform is defined as

.. math::
        f(r)=-\frac{1}{\pi}\int_r^\infty\frac{F'(y)}{\sqrt{y^2-r^2}}dy\approx-\frac{1}{\pi}\int_r^R\frac{F'(y)}{\sqrt{y^2-r^2}}dy\; .

**openAbel** takes care of taking the derivate of the input data supplied by the user.
The Backward Abel Transform is chosen  by setting :code:`forwardBackward=1`.



Backward (or Inverse) Abel Transform with Derivative Input
--------------
The Backward (or Inverse) Abel Transform with derivative input is defined as

.. math::
        f(r)=-\frac{1}{\pi}\int_r^\infty\frac{g(y)}{\sqrt{y^2-r^2}}dy\approx-\frac{1}{\pi}\int_r^R\frac{g(y)}{\sqrt{y^2-r^2}}dy\; .

In contrast to the normal Backward Abel Transform, **openAbel** expects to get the derivative as input by the user.

The Backward Abel Transform with derivative input is chosen by setting :code:`forwardBackward=2`.



Modified Forward Abel Transform
--------------

What we call the Modified Forward Abel Transform in **openAbel** is defined as the integral

.. math::
        H(y)=2\int_y^\infty\frac{h(r)y^2}{r^2\sqrt{r^2-y^2}}dr\approx2\int_y^R\frac{h(r)y^2}{r^2\sqrt{r^2-y^2}}dr\; .

I encountered this integral when a radial electric field of an atom (which has a :math:`1/r^2` singularity we want to integrate properly) is projected instead of a simpler function like with the normal Forward Abel Transform. One could just use the parameter :code:`shift = 0.5` instead to avoid the singularity of the electric field, but if one incorporates the singularity in the actual integral the convergence is much better. I recommend writing similar methods if one encounters other types of singularities in the Abel Transform.

The Modified Forward Abel Transform is chosen by setting :code:`forwardBackward=-2`.



