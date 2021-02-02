from matplotlib import pyplot as plt
from math import log, exp, sqrt
from scipy.special import erfc, erf
import matplotlib.colors as colors
import numpy as np

plt.style.use('latex')


POWER_STEP = 1.1
L_0_RANGE=[1.0e32, 2.0e34]
SIGMA_L_RANGE=[0.001, 1]

DIM_L0 = int(log(L_0_RANGE[1] / L_0_RANGE[0]) / log(POWER_STEP))
DIM_SIGMA = 50

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=1/5.0

DRAW_EXTRA_CONTOURS = False
DRAW_PLOEG_POINT = True

paperPoint = [0.88e34, 0.62]
ploegPoint = [10**32.206, 0.70585]

PATH_TO_FILE = "C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models-position/data/log-normal/"

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
        enterLine.append(float(item))
    lumSeen.append(np.asarray(enterLine))

totalNum = np.transpose(np.stack(totalNum, axis=0))
numSeen = np.transpose(np.stack(numSeen, axis=0))
lumSeen = np.transpose(np.stack(lumSeen, axis=0))

# ========================= Display data =========================

xVals = [L_0_RANGE[0] * POWER_STEP**i for i in range(DIM_L0)]
yVals = [SIGMA_L_RANGE[0] + (SIGMA_L_RANGE[1]-SIGMA_L_RANGE[0]) / DIM_SIGMA * j for j in range(DIM_SIGMA)]


fig, ax = plt.subplots(figsize=(7, 5))
plt.xlim(left=L_0_RANGE[0], right=L_0_RANGE[1])
plt.ylim(bottom=SIGMA_L_RANGE[0], top=SIGMA_L_RANGE[1])

plt.xscale("log")
plt.ylabel("$\sigma$")
plt.xlabel("$L_0$")
plt.title("Log normal luminosity function")

cols = colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                   vmax=max([max(v) for v in totalNum]))
c1 = plt.pcolor(xVals, yVals, totalNum, 
                   norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                   vmax=max([max(v) for v in totalNum])), cmap='Greys_r')
plt.colorbar(c1, extend='max')

# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(xVals, yVals, numSeen, [10*i for i in range(1, 10)], 
        colors=[(0, i/10.0, 0, 1) for i in range(1, 10)], linewidths=1)
plt.contour(xVals, yVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[(0, 0, 0)], linewidths=2, label="Number constraint")

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(xVals, yVals, lumSeen, [0.5 * i for i in range(1, 10)], 
        colors=[(1, i/10.0, 1-i/10.0, 1) for i in range(1, 10)], linewidths=1)
plt.contour(xVals, yVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[(0, 0, 0)], linestyles='dashed', linewidths=2, label="Fraction constraint")


# Plot thresholds
plt.plot(L_0_RANGE, [0.62-0.16, 0.62-0.16], c='blue', linewidth=1)
plt.plot(L_0_RANGE, [0.62+0.15, 0.62+0.15], c='blue', linewidth=1)
plt.plot([(0.88-0.41) * 1e34, (0.88-0.41) * 1e34], SIGMA_L_RANGE, c='blue', linewidth=1)
plt.plot([(0.88+0.79) * 1e34, (0.88+0.79) * 1e34], SIGMA_L_RANGE, c='blue', linewidth=1)

plt.scatter(paperPoint[0], paperPoint[1], c='blue')
#plt.scatter(minPoint[0], minPoint[1], c='cyan')

from matplotlib.lines import Line2D

# Ploeg point
if DRAW_PLOEG_POINT:
    plt.scatter(ploegPoint[0], ploegPoint[1], c='green')

custom_lines = [Line2D([0], [0], color='black', lw=2),
                Line2D([0], [0], color='black', lw=2, dashes=(4, 2))]
plt.legend(custom_lines, ['Number constraint', 'fraction constraint'])

plt.tight_layout()

# Save
if(DRAW_EXTRA_CONTOURS):
    plt.savefig(PATH_TO_FILE+"overlay-extra.png")
if(not DRAW_EXTRA_CONTOURS):
    plt.savefig(PATH_TO_FILE+"overlay.png")

plt.show()