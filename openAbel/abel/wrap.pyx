

import sys
import numpy

cimport base



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
        - '-2' modified Abel transform.
    shift : double
        Shift of the first sample away from 0 in units of stepSize.
        Usually this is either 0 or 0.5, and some methods only support these
        two values.
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
                 int method = 3, int order = 2, double eps = 1.e-15):
        
        try:
            self.plan = base.plan_fat(nData, forwardBackward, shift, stepSize, 
                                      method = method, order = order, eps = eps)
        except:
            print "Unexpected error in Cython routines:", sys.exc_info()[0], sys.exc_info()[1]
            raise


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
            '1' and '2' are not supported.

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

        dataInTemp = numpy.copy(dataIn)
        dataOut = numpy.copy(dataInTemp)

        try:
            base.execute_fat(self.plan, &dataInTemp[0], &dataOut[0], leftBoundary = leftBoundary, rightBoundary = rightBoundary)
        except:
            print "Unexpected error in Cython routines:", sys.exc_info()[0], sys.exc_info()[1]
            raise        

        return numpy.asarray(dataOut)[:self.plan.nData]


    def __dealloc__(self):

        base.destroy_fat(self.plan)

