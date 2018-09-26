
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
xMax = 1.
sig = 1./4.
stepSize = xMax/(nData-1)
noiseAmp = 0.01

abelObj = openAbel.Abel(nData, 2, 0., stepSize)  # Backward Abel transform where user inputs derivative


############################################################################################################################################
# No filtering
der = np.asarray([0.5, 0., -0.5])/stepSize
xx = np.linspace(-stepSize*(der.shape[0]-1)/2, xMax+stepSize*(der.shape[0]-1)/2, nData+(der.shape[0]-1))
dataIn = np.exp(-0.5*xx**2/sig**2)
np.random.seed(2202)
dataInWithNoise = dataIn + noiseAmp*np.random.randn(nData+(der.shape[0]-1))
dataOutAna = dataIn/np.sqrt(2*np.pi)/sig*erf(np.sqrt((xMax**2-xx**2)/2)/sig)

# Take derivatives
dataInD = np.convolve(dataInWithNoise, der, mode = 'valid')

# Backward transform
dataOutNoFilter = abelObj.execute(dataInD)


############################################################################################################################################
# Maximally flat filtering
# The length of this filter should be adjusted to the noise
der = np.asarray([4.76837e-7, 9.53674e-6, 0.0000901222, 0.000534058, 0.00221968, 0.00684929, 
                  0.0161719, 0.0295715, 0.041585, 0.0431252, 0.0280313, 0., -0.0280313, 
                  -0.0431252, -0.041585, -0.0295715, -0.0161719, -0.00684929,
                  -0.00221968, -0.000534058, -0.0000901222, -9.53674e-6, -4.76837e-7])/stepSize
xx = np.linspace(-stepSize*(der.shape[0]-1)/2, xMax+stepSize*(der.shape[0]-1)/2, nData+(der.shape[0]-1))
dataIn = np.exp(-0.5*xx**2/sig**2)
np.random.seed(2202)
dataInWithNoise = dataIn + noiseAmp*np.random.randn(nData+(der.shape[0]-1))
dataOutAna = dataIn/np.sqrt(2*np.pi)/sig*erf(np.sqrt((xMax**2-xx**2)/2)/sig)

# Take derivatives
dataIn = np.convolve(dataInWithNoise, der, mode = 'valid')

# Backward transform
dataOutMaxFlat = abelObj.execute(dataIn)

############################################################################################################################################
# Plotting and analytical (no noise) solution
xx = np.linspace(0., xMax, nData)
dataOutAna = np.exp(-0.5*xx**2/sig**2)/np.sqrt(2*np.pi)/sig*erf(np.sqrt((xMax**2-xx**2)/2)/sig)

fig, axarr = mpl.subplots(2, 1, sharex=True)

axarr[0].plot(xx, dataOutAna, 'r--', label='analytical')
axarr[0].plot(xx, dataOutNoFilter, 'b:', label='no filter')
axarr[0].plot(xx, dataOutMaxFlat, 'g-.', label='maxFlat')
axarr[0].set_ylabel('value')
axarr[0].legend()

axarr[1].semilogy(xx[:-1], np.abs((dataOutNoFilter[:-1]-dataOutAna[:-1])/dataOutAna[:-1]), 'b:', label='no filter')
axarr[1].semilogy(xx[:-1], np.abs((dataOutMaxFlat[:-1]-dataOutAna[:-1])/dataOutAna[:-1]), 'g-.', label='maxFlat')
axarr[1].set_ylabel('relative error')
axarr[1].set_xlabel('y')
axarr[1].legend()

fig.suptitle('Backward Abel Transform of a Gaussian', fontsize=16)

mpl.tight_layout()
mpl.subplots_adjust(top=0.9)



mpl.show()




