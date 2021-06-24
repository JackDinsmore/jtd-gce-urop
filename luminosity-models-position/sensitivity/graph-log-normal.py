from matplotlib import pyplot as plt
from math import log, exp, sqrt
from scipy.special import erfc, erf
import matplotlib.colors as colors
import numpy as np

from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import ImageGrid

plt.style.use('jcap')


PLOT_SIZE = 50
L_0_RANGE=[1.0e30, 2.0e36]#[1.0e32, 2.0e34]
SIGMA_L_RANGE=[0.001, 1]

lOPowerStep = (L_0_RANGE[1] / L_0_RANGE[0]) ** (1.0 / PLOT_SIZE)

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/5.0
TOTAL_FLUX = 1.2953417255755896e-09#7.494712733226778e-10

DRAW_EXTRA_CONTOURS = False
LINE_COLOR = "C1"
DRAW_PLOEG_POINT = True

paperPoint = [0.88e34, 0.62]
ploegPoint = [10**32.206, 0.70585]
gautamPoint = [3.91983577e+32, 0.937184991]
SHADE_SCALE=25

def get_path(mult):
    return "/home/jtdinsmo/Dropbox (MIT)/GCE UROP/luminosity-models-position/data-"\
        + str(mult) + "x/log-normal/"

def shade(ax, field, threshold, xs, ys, off=False):
    px = []
    py = []
    for x in range(0 if off else 1, SHADE_SCALE, 1):
        inx = int(float(x) / SHADE_SCALE * field.shape[1])
        for y in range(0 if off else 1, SHADE_SCALE, 1):
            iny = int(float(y) / SHADE_SCALE * field.shape[0])
            if field[iny][inx] < threshold:
                fracx = float(x) / SHADE_SCALE * field.shape[1] - inx
                fracy = float(y) / SHADE_SCALE * field.shape[0] - iny
                px.append(xs[inx] + fracx * (xs[inx+1] - xs[inx]))
                py.append(ys[iny] + fracy * (ys[iny+1] - ys[iny]))
    ax.scatter(px, py, marker=('|' if off else '_'), c=LINE_COLOR, sizes = (20,), alpha=0.7)

# ========================== Load data ===========================


xVals = [L_0_RANGE[0] * lOPowerStep**i for i in range(PLOT_SIZE)]
yVals = [SIGMA_L_RANGE[0] + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / PLOT_SIZE * j for j in range(PLOT_SIZE)]

fig = plt.figure(figsize=(4.5, 12))
axs = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(4,1),
                 axes_pad=0.15,
                 aspect=False,
                 share_all=False,
                 cbar_location="bottom",
                 cbar_mode="single",
                 cbar_size="7%",
                 cbar_pad=0.6,
                 )

for ax in axs:
    ax.set_xlim(left=min(xVals), right=max(xVals))
    ax.set_ylim(bottom=min(yVals), top=max(yVals))
    ax.set_ylabel("$\sigma$")
    ax.set_xscale("log")
axs[-1].set_xlabel("$L_0$ [erg / s]")

i=0
for mult in [1, 2, 5, 10]:
    totalNum = []
    numSeen = []
    lumSeen = []

    f = open(get_path(mult) + "total-num.txt")
    for line in f.read().split('\n')[:-1]:
        enterLine = []
        for item in line.split(', '):
            if item == '' or item == ' ': continue
            enterLine.append(float(item))
        totalNum.append(np.asarray(enterLine))

    f = open(get_path(mult) + "num-seen.txt")
    for line in f.read().split('\n')[:-1]:
        enterLine = []
        for item in line.split(', '):
            if item == '' or item == ' ': continue
            enterLine.append(float(item))
        numSeen.append(np.asarray(enterLine))

    f = open(get_path(mult) + "lum-seen.txt")
    for line in f.read().split('\n')[:-1]:
        enterLine = []
        for item in line.split(', '):
            if item == '' or item == ' ': continue
            enterLine.append(float(item) / TOTAL_FLUX)
        lumSeen.append(np.asarray(enterLine))

    totalNum = np.transpose(np.stack(totalNum, axis=0))
    numSeen = np.transpose(np.stack(numSeen, axis=0))
    lumSeen = np.transpose(np.stack(lumSeen, axis=0))

    # ========================= Display data =========================
    c1 = axs[i].contourf(xVals, yVals, totalNum,
        norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
            vmax=max([max(v) for v in totalNum])), cmap='Greys_r')

    axs[i].contour(xVals, yVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[LINE_COLOR])
    axs[i].contour(xVals, yVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed')

    axs[i].plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='^', markersize=6)
    axs[i].errorbar([paperPoint[0]], [paperPoint[1]], xerr=[[0.41e34], [0.79e34]], yerr=[[0.16], [0.15]], linewidth=1, color=LINE_COLOR)
    axs[i].plot(ploegPoint[0], ploegPoint[1], markeredgecolor='black', markerfacecolor="fuchsia", marker='s', markersize=6)
    axs[i].errorbar([paperPoint[0]], [paperPoint[1]], xerr=[[4.1304e+31], [4.1304e+31]], yerr=[[0.01011908], [0.01011908]], linewidth=1, color="fuchsia")
    axs[i].plot(gautamPoint[0], gautamPoint[1], markeredgecolor='black', markerfacecolor="red", marker='*', markersize=8)
    axs[i].errorbar([paperPoint[0]], [paperPoint[1]], xerr=[[3.0842e+31], [3.0842e+31]], yerr=[[0.00316971], [0.00316971]], linewidth=1, color="red")

    # Observation
    shade(axs[i], numSeen, NUM_PULSARS_ABOVE_THRESHOLD, xVals, yVals)
    shade(axs[i], lumSeen, FRAC_ABOVE_THRESHOLD, xVals, yVals, True)
    i+=1


cbar = axs[-1].cax.colorbar(c1)
cbar.set_label("$N_\\textrm{GCE}$")

custom_lines = [Line2D([0], [0], color=LINE_COLOR),
                Line2D([0], [0], color=LINE_COLOR, dashes=(4, 2)),
                Line2D([], [], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='^', linestyle='None', markersize=6),
                Line2D([], [], markeredgecolor='black', markerfacecolor="C6", marker='s', linestyle='None', markersize=6),
                Line2D([], [], markeredgecolor='black', markerfacecolor="red", marker='*', linestyle='None', markersize=8),]
axs[-1].legend(custom_lines,
    ['$N_\\textrm{r} = 47$', '$R_\\textrm{r}=0.2$', "GLC", "GCE", "AIC"],
    loc="lower left")
fig.tight_layout()

plt.tight_layout()
plt.savefig("log-normal-sensitivity.pdf")
plt.show()
