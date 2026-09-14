############################################################################################################################################
# Simple example which shows how to select different methods and orders.
# Results are compared with the analytical solution.
############################################################################################################################################


import matplotlib.pyplot as mpl
import numpy as np
from scipy.special import erf

import openabel

############################################################################################################################################
# Plotting setup
# This block can be ignored, it's just for nicer plots.

params = {
    "axes.labelsize": 8,
    "font.size": 8,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": False,
    "figure.figsize": [5.0, 5.0],
}
mpl.rcParams.update(params)
# Color scheme
colors = [
    "#005AA9",
    "#E6001A",
    "#99C000",
    "#721085",
    "#EC6500",
    "#009D81",
    "#A60084",
    "#0083CC",
    "#F5A300",
    "#C9D400",
    "#FDCA00",
]
# Plot markers
markers = ["o", "v", "s", "D", "p", "*", "h", "+", "^", "x"]
# Line styles
linestyles = ["-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":"]
lw = 2

############################################################################################################################################

# Parameters
n_data = 40
shift = 0.0
x_max = 3.5
sig = 1.0
step_size = x_max / (n_data - 1)

forward_backward = -1  # Forward transform, similar definition ('1' = backward) as in FFT libraries.

# Create Abel transform object for three different methods and orders.
# Some methods ignore the order keyword argument, and for the normal user
# only method = 3 and order = 2 to order = 5 are recommended.
# Higher orders require data outside the integration domain to be stable.
# For more information see the documentation.
abel_obj0 = openabel.Abel(n_data, forward_backward, shift, step_size, method=2, order=2)
abel_obj1 = openabel.Abel(n_data, forward_backward, shift, step_size, method=3, order=5)
abel_obj2 = openabel.Abel(n_data, forward_backward, shift, step_size, method=3, order=11)

# Input data
xx = np.linspace(shift * step_size, x_max, n_data)
data_in = np.exp(-0.5 * xx**2 / sig**2)
xx_ext = np.linspace(
    shift * step_size,
    x_max + 5 * step_size,
    n_data + 5,
)  # floor((order-1)/2) extra points at right end
data_in_ext = np.exp(-0.5 * xx_ext**2 / sig**2)

# Backward transform and analytical result

data_out0 = abel_obj0.execute(data_in)
data_out1 = abel_obj1.execute(data_in)
data_out2 = abel_obj2.execute(
    data_in_ext,
    left_boundary=2,
    right_boundary=3,
)  # 2 means use even symmetry, 3 means input extra points.
data_out_ana = data_in * np.sqrt(2 * np.pi) * sig * erf(np.sqrt((x_max**2 - xx**2) / 2) / sig)

# Plotting
fig, axarr = mpl.subplots(2, 1, sharex=True)

axarr[0].plot(xx, data_out_ana, color=colors[0], marker=markers[0], linestyle=linestyles[0], label="analy.")
axarr[0].plot(xx, data_out0, color=colors[1], marker=markers[1], linestyle=linestyles[1], label="openAbel TE 2nd")
axarr[0].plot(xx, data_out1, color=colors[2], marker=markers[2], linestyle=linestyles[2], label="openAbel FMM 5th")
axarr[0].plot(xx, data_out2, color=colors[3], marker=markers[3], linestyle=linestyles[3], label="openAbel FMM 11th")
axarr[0].set_ylabel("value")
axarr[0].legend()

axarr[1].semilogy(
    xx[:-1] / sig,
    np.abs((data_out0[:-1] - data_out_ana[:-1]) / data_out_ana[:-1]),
    color=colors[1],
    marker=markers[1],
    linestyle=linestyles[1],
    label="openAbel TE 2nd",
)
axarr[1].semilogy(
    xx[:-1] / sig,
    np.abs((data_out1[:-1] - data_out_ana[:-1]) / data_out_ana[:-1]),
    color=colors[2],
    marker=markers[2],
    linestyle=linestyles[2],
    label="openAbel FMM 5th",
)
axarr[1].semilogy(
    xx[:-1] / sig,
    np.abs((data_out2[:-1] - data_out_ana[:-1]) / data_out_ana[:-1]),
    color=colors[3],
    marker=markers[3],
    linestyle=linestyles[3],
    label="openAbel FMM 11th",
)
axarr[1].set_ylabel("relative error")
axarr[1].set_xlabel("y")
axarr[1].legend()

mpl.tight_layout()
mpl.savefig("example002_method_order.png", dpi=300)

mpl.show()
