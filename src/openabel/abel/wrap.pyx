

import numpy

cimport openabel.abel.base as base
cimport openabel.constants as const

cdef class Abel(object):
    """
    This is a wrapper class to make the Abel transform from Cython available
    in Python.

    Parameters
    ----------
    n_data : int
        Length of the data vector.
    forward_backward : int
        Which transform to perform:
        - '-1' forward Abel transform
        - '1' backward (or inverse) Abel transform
        - '2' backward (or inverse) Abel transform with the 
          derivative already supplied by user
        - '-2' modified forward Abel transform.
    shift : double
        Shift of the first sample away from 0 in positive direction in units of step_size.
        Usually this is either 0 or 0.5, and some methods only support these two values.
    step_size : double
        Step size (or grid spacing) between two data points.
    method : int, optional
        Which method to employ to calculate transform.

    Raises
    ------
    ValueError
        If an input parameter has a not viable value.
    NotImplementedError
        If a method doesn't (yet) support the operation given by parameters.
    """

    def __init__(self, int n_data, int forward_backward, double shift, double step_size, 
                 int method = 3, int order = 2, double eps = 1.e1*const.machine_epsilon):

        cdef int order_filter

        if n_data < 2:
            raise ValueError('n_data must be at least 2.')

        try:
            self.plan = base.plan_fat(n_data, forward_backward, shift, step_size, 
                                      method = method, order = order, eps = eps)
        except:
            raise

        # Samples per side that execute() reads beyond n_data with boundary value 3: the half widths of the
        # end-correction stencil and of the derivative filter, the same extension widths trap.pyx and fmm.pyx use.
        # Method 0 has a first-order stencil, method 1 ignores the boundary values.
        if method == 1:
            self.n_outside = 0
        else:
            if method == 0:
                order = 1
            order_filter = order+1 + (order % 2) if forward_backward == 1 else 1
            self.n_outside = (order-1)//2 + (order_filter-1)//2


    # TODO maybe support 2D (or nD) arrays as well here?
    def execute(self, double[:] data_in, int left_boundary = 0, int right_boundary = 0):
        """
        This is the function which actually does the transform.

        Parameters
        ----------
        data_in : numpy.array
            Data vector.
        left_boundary : int, optional
            Defines how the start of the data are handled:
            - '0' data only given inside integration interval
            - '1' data has odd symmetry around zero
            - '2' data has even symmetry around zero
            - '3' data is given outside domain as well.
        right_boundary : int, optional
            Almost the same as `left_boundary` only for end the data.
            '1' and '2' are not supported here.

        Returns
        ------
        data_out : numpy.array
            Transformed data.
            
        Raises
        ------
        ValueError
            If an input parameter has a not viable value.
        NotImplementedError
            If a method doesn't (yet) support the operation given by parameters.
        """
        cdef:
            double[::1] data_in_temp
            double[::1] data_out
            Py_ssize_t n_needed = self.plan.n_data

        if left_boundary == 3:
            n_needed += self.n_outside
        if right_boundary == 3:
            n_needed += self.n_outside
        if data_in.shape[0] < n_needed:
            raise ValueError(f'data_in has {data_in.shape[0]} samples, but the plan needs at least {n_needed} for the '
                             f'given boundary values.')

        data_in_temp = numpy.copy(data_in)
        data_out = numpy.copy(data_in_temp)

        try:
            base.execute_fat(self.plan, &data_in_temp[0], &data_out[0], left_boundary = left_boundary, 
                             right_boundary = right_boundary)
        except:
            raise        

        return numpy.asarray(data_out)[:self.plan.n_data]


    def __dealloc__(self):

        base.destroy_fat(self.plan)
