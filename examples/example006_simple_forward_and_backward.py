############################################################################################################################################
# Simple example which calculates forward and backward Abel transform of a Gaussian.
# Results are compared with the analytical solution. Mostly default parameters are used.
############################################################################################################################################


import matplotlib.pyplot as mpl
import numpy as np

import openabel

############################################################################################################################################
# Plotting setup
# This block can be ignored, it's just for nicer plots.

params = {
    "axes.labelsize": 8,
    "font.size": 8,
    "legend.fontsize": 8,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": False,
    "figure.figsize": [6.5, 5.0],
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
n_data = 100000
shift = 0.0
x_max = 20.0
sig = 1.0
step_size = x_max / (n_data - 1)

# Create Abel transform object, which does all precomputation possible without knowing the exact data.
abel_obj_fw = openabel.Abel(n_data, -1, shift, step_size, order=3)
abel_obj_bw = openabel.Abel(n_data, -1, shift, step_size, order=3)

# Input data
xx = np.linspace(shift * step_size, x_max, n_data)
data_in = np.exp(-0.5 * xx**2 / sig**2)

# Forward Abel transform and analytical result.
# We show both the analytical result of a truncated Gaussian and a standard Gaussian to show
# that some error is due to truncation.
data_out = abel_obj_fw.execute(data_in)
data_out = abel_obj_fw.execute(data_out)
for ii in range(n_data):
    data_out[ii] /= 2.0 * np.pi


# Plotting
fig, axarr = mpl.subplots(2, 1, sharex=True, layout="constrained")

axarr[0].plot(
    xx,
    data_in,
    color=colors[0],
    marker=markers[0],
    linestyle=linestyles[0],
    markevery=n_data // 40,
    label="analy.",
)
axarr[0].plot(
    xx,
    data_out,
    color=colors[2],
    marker=markers[2],
    linestyle=linestyles[2],
    markevery=n_data // 40,
    label="openAbel",
)
axarr[0].set_ylabel("value")
axarr[0].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0))

axarr[1].semilogy(
    xx[:-1],
    np.abs((data_out[:-1] - data_in[:-1]) / data_in[:-1]),
    color=colors[3],
    marker=markers[3],
    linestyle=linestyles[3],
    markevery=n_data // 40,
    label="not trunc.",
)
axarr[1].set_ylabel("relative error")
axarr[1].set_xlabel("y")
axarr[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0))

mpl.savefig("example006_simple_forward_and_backward.png", dpi=300)

mpl.show()
