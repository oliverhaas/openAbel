

import numpy

cimport openabel.abel.base as base
cimport openabel.constants as const

cdef class Abel(object):
    """
    This is a wrapper class to make the Abel transform from Cython available
    in Python.

    Parameters
    ----------
    nData : int
        Length of the data vector.
    forwardBackward : int
        Which transform to perform:
        - '-1' forward Abel transform
        - '1' backward (or inverse) Abel transform
        - '2' backward (or inverse) Abel transform with the 
          derivative already supplied by user
        - '-2' modified forward Abel transform.
    shift : double
        Shift of the first sample away from 0 in positive direction in units of stepSize.
        Usually this is either 0 or 0.5, and some methods only support these two values.
    stepSize : double
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

    def __init__(self, int nData, int forwardBackward, double shift, double stepSize, 
                 int method = 3, int order = 2, double eps = 1.e1*const.machineEpsilon):

        cdef int orderFilter

        if nData < 2:
            raise ValueError('nData must be at least 2.')

        try:
            self.plan = base.plan_fat(nData, forwardBackward, shift, stepSize, 
                                      method = method, order = order, eps = eps)
        except:
            raise

        # Samples per side that execute() reads beyond nData with boundary value 3: the half widths of the
        # end-correction stencil and of the derivative filter, the same extension widths trap.pyx and fmm.pyx use.
        # Method 0 has a first-order stencil, method 1 ignores the boundary values.
        if method == 1:
            self.nOutside = 0
        else:
            if method == 0:
                order = 1
            orderFilter = order+1 + (order % 2) if forwardBackward == 1 else 1
            self.nOutside = (order-1)//2 + (orderFilter-1)//2


    # TODO maybe support 2D (or nD) arrays as well here?
    def execute(self, double[:] dataIn, int leftBoundary = 0, int rightBoundary = 0):
        """
        This is the function which actually does the transform.

        Parameters
        ----------
        dataIn : numpy.array
            Data vector.
        leftBoundary : int, optional
            Defines how the start of the data are handled:
            - '0' data only given inside integration interval
            - '1' data has odd symmetry around zero
            - '2' data has even symmetry around zero
            - '3' data is given outside domain as well.
        rightBoundary : int, optional
            Almost the same as `leftBoundary` only for end the data.
            '1' and '2' are not supported here.

        Returns
        ------
        dataOut : numpy.array
            Transformed data.
            
        Raises
        ------
        ValueError
            If an input parameter has a not viable value.
        NotImplementedError
            If a method doesn't (yet) support the operation given by parameters.
        """
        cdef:
            double[::1] dataInTemp
            double[::1] dataOut
            Py_ssize_t nNeeded = self.plan.nData

        if leftBoundary == 3:
            nNeeded += self.nOutside
        if rightBoundary == 3:
            nNeeded += self.nOutside
        if dataIn.shape[0] < nNeeded:
            raise ValueError(f'dataIn has {dataIn.shape[0]} samples, but the plan needs at least {nNeeded} for the '
                             f'given boundary values.')

        dataInTemp = numpy.copy(dataIn)
        dataOut = numpy.copy(dataInTemp)

        try:
            base.execute_fat(self.plan, &dataInTemp[0], &dataOut[0], leftBoundary = leftBoundary, 
                             rightBoundary = rightBoundary)
        except:
            raise        

        return numpy.asarray(dataOut)[:self.plan.nData]


    def __dealloc__(self):

        base.destroy_fat(self.plan)
