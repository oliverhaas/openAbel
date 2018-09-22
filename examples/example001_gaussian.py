############################################################################################################################################
# Simple example which calculates forward and backward Abel transform of a Gaussian.
# Results are compared with the analytical solution. Mostly default parameters are used.
############################################################################################################################################


import openAbel
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as mpl


# Parameters
nData = 200
forwardBackward = -1    # Forward transform, similar definition to many FFT libraries.
shift = 0.
xMax = 4.2
sig = 1.
stepSize = xMax/(nData-1)

# Create forward Abel transform object, which does all precomputation possible without knowing the exact 
# data. This way it's much faster if repeated transforms are done.
abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize)    # Default keyword parameters are (method = 4, order = 2, eps = 1.e-15)

# Input data
xx = np.linspace(shift*stepSize, xMax, nData)
dataIn = np.exp(-0.5*xx**2/sig**2)

# Forward Abel transform and analytical result
# We show both the analytical result of a truncated Gaussian and a standard Gaussian to show
# that some error is due to truncation.
dataOut = abelObj.execute(dataIn)
dataOutAna = dataIn*np.sqrt(2*np.pi)*sig
dataOutAnaTrunc = dataIn*np.sqrt(2*np.pi)*sig*erf(np.sqrt((xMax**2-xx**2)/2)/sig)


# Plotting
fig, axarr = mpl.subplots(2, 1, sharex=True)

axarr[0].plot(xx, dataOutAna, 'r--', label='analy.')
axarr[0].plot(xx, dataOutAnaTrunc, 'g-.', label='analy. trunc.')
axarr[0].plot(xx, dataOut, 'b:', label='openAbel')
axarr[0].set_ylabel('value')
axarr[0].legend()

axarr[1].semilogy(xx[:-1]/sig, np.abs((dataOut[:-1]-dataOutAna[:-1])/dataOutAna[:-1]), 'r--', label='rel. err.')
axarr[1].semilogy(xx[:-1]/sig, np.abs((dataOut[:-1]-dataOutAnaTrunc[:-1])/dataOutAnaTrunc[:-1]), 'b:', label='rel. err. trunc.')
axarr[1].set_ylabel('relative error')
axarr[1].set_xlabel('y')
axarr[1].legend()

fig.suptitle('Forward Abel Transform of a Gaussian', fontsize=16)

mpl.tight_layout()
mpl.subplots_adjust(top=0.87)
mpl.show(block = False)



############################################################################################################################################
# Now the same thing for the inverse transform
forwardBackward = 1
abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize) # Default keyword parameters are (method = 4, order = 2, eps = 1.e-15)

# Backward transform
dataOut = abelObj.execute(dataIn)
dataOutAna = dataIn/np.sqrt(2*np.pi)/sig
dataOutAnaTrunc = dataIn/np.sqrt(2*np.pi)/sig*erf(np.sqrt((xMax**2-xx**2)/2)/sig)

# Plotting
fig, axarr = mpl.subplots(2, 1)

axarr[0].plot(xx, dataOutAna, 'r--', label='analy.')
axarr[0].plot(xx, dataOutAnaTrunc, 'g-.', label='analy. trunc.')
axarr[0].plot(xx, dataOut, 'b:', label='openAbel')
axarr[0].set_title('Value')
axarr[0].legend()

axarr[1].semilogy(xx[:-1], np.abs((dataOut[:-1]-dataOutAna[:-1])/dataOutAna[:-1]), 'r--', label='rel. err.')
axarr[1].semilogy(xx[:-1], np.abs((dataOut[:-1]-dataOutAnaTrunc[:-1])/dataOutAnaTrunc[:-1]), 'b:', label='rel. err. trunc.')
axarr[1].set_title('Relative error')
axarr[1].legend()

fig.suptitle('Backward Abel Transform of a Gaussian', fontsize=16)

mpl.tight_layout()
mpl.subplots_adjust(top=0.87)

#fig, ax = mpl.subplots(1, 1)
#for order in range(1,4):
#    for method in range(3,5):

#        forwardBackward = 1
#        abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize, method = method, order = order)

#        dataOut = abelObj.execute(dataIn)
#        dataOutAna = dataIn/np.sqrt(2*np.pi)/sig
#        dataOutAnaTrunc = dataIn/np.sqrt(2*np.pi)/sig*erf(np.sqrt((xMax**2-xx**2)/2)/sig)

#        ax.semilogy(xx[:-1], np.abs((dataOut[:-1]-dataOutAnaTrunc[:-1])/dataOutAnaTrunc[:-1]), label=str(method)+','+str(order))
#        ax.legend()

#        mpl.tight_layout()
#        mpl.subplots_adjust(top=0.87)

mpl.show()
