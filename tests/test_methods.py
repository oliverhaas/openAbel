

import openAbel as oa
from nose.tools import *
import numpy as np
from scipy.special import erf
import itertools as it


@nottest
def helper_test_method(nData, forwardBackward, shift, methods, orders, rtol):

    for method, order in it.product(methods,orders):
        xMax = 3.5
        sig = 1.
        stepSize = xMax/(nData-1)

        abelObj = oa.Abel(nData, forwardBackward, shift, stepSize, method = method, order = order)    

        xx = np.linspace(shift*stepSize, xMax, nData)

        dataIn = np.exp(-0.5*xx**2/sig**2)
        if forwardBackward == -1:
            dataOutAna = dataIn*np.sqrt(2*np.pi)*sig*erf(np.sqrt((xMax**2-xx**2)/2)/sig)
        elif forwardBackward == 1:
            dataOutAna = dataIn/np.sqrt(2*np.pi)/sig*erf(np.sqrt((xMax**2-xx**2)/2)/sig)
        elif forwardBackward == 2:
            dataOutAna = dataIn/np.sqrt(2*np.pi)/sig*erf(np.sqrt((xMax**2-xx**2)/2)/sig)
            dataIn = -xx/sig**2*np.exp(-0.5*xx**2/sig**2)
        else:
            raise NotImplementedError('Test not implemented.')
        
        dataOut = abelObj.execute(dataIn) 
        
        assert dataOut[-1] == 0.
        np.testing.assert_allclose(dataOut[:-1],dataOutAna[:-1],rtol=rtol)

    return None


def test_methodEndCorr():
    helper_test_method(200, -1, 0., [2,3], np.arange(1,6), 1.e-2)
    helper_test_method(200, 1, 0., [2,3], np.arange(1,6), 1.e-1)
    helper_test_method(200, 2, 0., [2,3], np.arange(1,6), 1.e-2)
    

def test_hansenLaw():
    helper_test_method(200, -1, 0., [1], [1], 1.e-1)
    
    
def test_trapDesing():
    helper_test_method(200, -1, 0., [0], [1], 1.e-1)
