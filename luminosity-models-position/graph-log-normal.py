from matplotlib import pyplot as plt
from math import log, exp, sqrt
from scipy.special import erfc, erf
import matplotlib.colors as colors
import numpy as np

from matplotlib.lines import Line2D

plt.style.use('jcap')


PLOT_SIZE = 50
L_0_RANGE=[1.0e30, 2.0e36]#[1.0e32, 2.0e34]
SIGMA_L_RANGE=[0.001, 1]

MULTIPLIER = 1

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

PATH_TO_FILE = "/home/jtdinsmo/Dropbox (MIT)/GCE UROP/luminosity-models-position/data-"\
 + str(MULTIPLIER) + "x/log-normal/"
SHADE_SCALE=25

def shade(field, threshold, xs, ys, off=False):
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
    plt.scatter(px, py, marker=('|' if off else '_'), c=LINE_COLOR, sizes = (20,), alpha=0.7)

# ========================== Load data ===========================

totalNum = []
numSeen = []
lumSeen = []

f = open(PATH_TO_FILE + "total-num.txt")
for line in f.read().split('\n')[:-1]:
    enterLine = []
    for item in line.split(', '):
        if item == '' or item == ' ': continue
        enterLine.append(float(item))
    totalNum.append(np.asarray(enterLine))

f = open(PATH_TO_FILE + "num-seen.txt")
for line in f.read().split('\n')[:-1]:
    enterLine = []
    for item in line.split(', '):
        if item == '' or item == ' ': continue
        enterLine.append(float(item))
    numSeen.append(np.asarray(enterLine))

f = open(PATH_TO_FILE + "lum-seen.txt")
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

xVals = [L_0_RANGE[0] * lOPowerStep**i for i in range(PLOT_SIZE)]
yVals = [SIGMA_L_RANGE[0] + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / PLOT_SIZE * j for j in range(PLOT_SIZE)]


fig, ax = plt.subplots()
plt.xlim(left=L_0_RANGE[0], right=L_0_RANGE[1])
plt.ylim(bottom=SIGMA_L_RANGE[0], top=SIGMA_L_RANGE[1])

plt.xscale("log")
plt.ylabel("$\sigma$")
plt.xlabel("$L_0$ [ergs / s]")

cols = colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                   vmax=max([max(v) for v in totalNum]))
c1 = plt.contourf(xVals, yVals, totalNum,
                   norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                   vmax=max([max(v) for v in totalNum])), cmap='Greys_r')

cbar = plt.colorbar(c1, extend='max')
cbar.set_label("$N_\\textrm{GCE}$")


# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(xVals, yVals, numSeen, [10*i for i in range(1, 10)],
        colors=[(0, i/10.0, 0, 1) for i in range(1, 10)])
plt.contour(xVals, yVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[LINE_COLOR])

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(xVals, yVals, lumSeen, [0.5 * i for i in range(1, 10)],
        colors=[(1, i/10.0, 1-i/10.0, 1) for i in range(1, 10)])
plt.contour(xVals, yVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed')

plt.plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='^', markersize=6)
plt.errorbar([paperPoint[0]], [paperPoint[1]], xerr=[[0.41e34], [0.79e34]], yerr=[[0.16], [0.15]], linewidth=1, color=LINE_COLOR)
plt.plot(ploegPoint[0], ploegPoint[1], markeredgecolor='black', markerfacecolor="fuchsia", marker='s', markersize=6)
plt.errorbar([paperPoint[0]], [paperPoint[1]], xerr=[[4.1304e+31], [4.1304e+31]], yerr=[[0.01011908], [0.01011908]], linewidth=1, color="fuchsia")
plt.plot(gautamPoint[0], gautamPoint[1], markeredgecolor='black', markerfacecolor="red", marker='*', markersize=8)
plt.errorbar([paperPoint[0]], [paperPoint[1]], xerr=[[3.0842e+31], [3.0842e+31]], yerr=[[0.00316971], [0.00316971]], linewidth=1, color="red")

# Observation
shade(numSeen, NUM_PULSARS_ABOVE_THRESHOLD, xVals, yVals)
shade(lumSeen, FRAC_ABOVE_THRESHOLD, xVals, yVals, True)


# Final points


custom_lines = [Line2D([0], [0], color=LINE_COLOR),
                Line2D([0], [0], color=LINE_COLOR, dashes=(4, 2)),
                Line2D([], [], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='^', linestyle='None', markersize=6),
                Line2D([], [], markeredgecolor='black', markerfacecolor="C6", marker='s', linestyle='None', markersize=6),
                Line2D([], [], markeredgecolor='black', markerfacecolor="red", marker='*', linestyle='None', markersize=8),]
plt.legend(custom_lines, ['$N_\\textrm{r} = 47$', '$R_\\textrm{r}=0.2$', "GCL", "GCE", "AIC"], loc="lower left")
plt.xlim(xVals[0], xVals[-1])
plt.ylim(yVals[0], yVals[-1])

plt.tight_layout()

# Save
if(DRAW_EXTRA_CONTOURS):
    plt.savefig(PATH_TO_FILE+"overlay-extra.png")
if(not DRAW_EXTRA_CONTOURS):
    plt.savefig("log-normal-pos-x" + str(MULTIPLIER) + ".pdf")

plt.show()
