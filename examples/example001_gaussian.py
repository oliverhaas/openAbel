###################################################################################################
# Simple example which calculates forward and backward abel transform of a Gaussian.
# Results are compared with the analytical solution. Mostly default parameters are used.
#
###################################################################################################



import openAbel
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as mpl


# Parameters
nData = 1000
forwardBackward = -1    # Forward transform, similar definition to many FFT libraries.
shift = 0.
xMax = 8.
sig = 1.
stepSize = xMax/(nData-1)

# Create Abel transform object, which does all precomputation possible without knowing the exact 
# data. This way it's much faster if repeated transforms are done.
abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize) #  ,method = 4, order = 1, eps = 1.e-15)

# Input data
xx = np.linspace(shift*stepSize, xMax, nData)
dataIn = np.exp(-0.5*xx**2/sig**2)

# Forward Abel transform and analytical result
# We show both the analytical result of a truncated Gaussian and a standard Gaussian to show
# that some error is due to truncation.
dataOut = abelObj.execute(dataIn)
dataOutAna = dataIn*np.sqrt(2*np.pi)*sig
dataOutAnaTrunc = dataIn*np.sqrt(2*np.pi)*sig*erf(np.sqrt((xMax**2-xx**2)/2)/sig)


fig, axarr = mpl.subplots(2, 1)

axarr[0].plot(xx, dataOutAna, 'r--')
axarr[0].plot(xx, dataOutAnaTrunc, 'g-.')
axarr[0].plot(xx, dataOut, 'b:')

axarr[1].semilogy(xx, np.abs((dataOut-dataOutAna)/dataOutAna), 'r--')
axarr[1].semilogy(xx, np.abs((dataOut-dataOutAnaTrunc)/dataOutAnaTrunc), 'b:')


fig.suptitle('Forward Abel Transform of a Gaussian', fontsize=16)

mpl.tight_layout()
mpl.subplots_adjust(top=0.9)
mpl.show()
