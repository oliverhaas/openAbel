############################################################################################################################################
# Simple example which calculates backward Abel transform of a Gaussian.
# Results are compared with the analytical solution. Mostly default parameters are used.
############################################################################################################################################


import openAbel
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as mpl

############################################################################################################################################
# Plotting setup
# This block can be ignored, it's just for nicer plots.

params = {
   'axes.labelsize': 8,
   'font.size': 8,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'figure.figsize': [5., 5.]
   }
mpl.rcParams.update(params)
# Color scheme
colors = ['#005AA9','#E6001A','#99C000','#721085','#EC6500','#009D81','#A60084','#0083CC','#F5A300','#C9D400','#FDCA00']
# Plot markers
markers = ["o", "v" , "s", "D", "p", "*", "h", "+", "^", "x"]
# Line styles 
linestyles = ['-', '--', '-.', ':','-', '--', '-.', ':','-', '--', '-.', ':']
lw = 2

############################################################################################################################################

# Parameters
nData = 40
shift = 0.
xMax = 3.5
sig = 1.
stepSize = xMax/(nData-1)
forwardBackward = 1    # Backward transform, similar definition ('1' = backward) as in FFT libraries.

# Create Abel transform object, which does all precomputation possible without knowing the exact data.
abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize)    

# Input data
xx = np.linspace(shift*stepSize, xMax, nData)
dataIn = np.exp(-0.5*xx**2/sig**2)

# Backward transform and analytical result.
# We show both the analytical result of a truncated Gaussian and a standard Gaussian to show
# that some error is due to truncation.
dataOut = abelObj.execute(dataIn)
dataOutAna = dataIn/np.sqrt(2*np.pi)/sig
dataOutAnaTrunc = dataIn/np.sqrt(2*np.pi)/sig*erf(np.sqrt((xMax**2-xx**2)/2)/sig)

# There is the option for the user to provide the derivative in the backward Abel transform directly.
# This is useful and can decrease the error, e.g. if the derivative can be taken analytically.
forwardBackward = 2
abelObj = openAbel.Abel(nData, forwardBackward, shift, stepSize) 
dataIn = -xx/sig**2*np.exp(-0.5*xx**2/sig**2)
dataOut2 = abelObj.execute(dataIn)


# Plotting
fig, axarr = mpl.subplots(2, 1, sharex=True)

axarr[0].plot(xx, dataOutAna, color = colors[0], marker = markers[0], linestyle = linestyles[0], label='analy.')
axarr[0].plot(xx, dataOutAnaTrunc, color = colors[1], marker = markers[1], linestyle = linestyles[1], label='analy. trunc.')
axarr[0].plot(xx, dataOut, color = colors[2], marker = markers[2], linestyle = linestyles[2], label='openAbel')
axarr[0].plot(xx, dataOut2, color = colors[3], marker = markers[3], linestyle = linestyles[3], label='openAbel analy. der.')
axarr[0].set_ylabel('value')
axarr[0].legend()

axarr[1].semilogy(xx[:-1], np.abs((dataOut[:-1]-dataOutAna[:-1])/dataOutAna[:-1]),
                  color = colors[4], marker = markers[4], linestyle = linestyles[4], label = 'not trunc.')
axarr[1].semilogy(xx[:-1], np.abs((dataOut[:-1]-dataOutAnaTrunc[:-1])/dataOutAnaTrunc[:-1]),
                  color = colors[5], marker = markers[5], linestyle = linestyles[5],  label='trunc.')
axarr[1].semilogy(xx[:-1], np.abs((dataOut2[:-1]-dataOutAnaTrunc[:-1])/dataOutAnaTrunc[:-1]),
                  color = colors[6], marker = markers[6], linestyle = linestyles[6],  label='trunc. analy. der.')
axarr[1].set_ylabel('relative error')
axarr[1].set_xlabel('y')
axarr[1].legend()

mpl.tight_layout()
mpl.savefig('example001_simpleBackward.png', dpi=300)

mpl.show()



