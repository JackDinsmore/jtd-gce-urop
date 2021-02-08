''' Things to check:
    Is the NFW profile a number density?
    What is r_s from the Fermilab paper?
    I should integrate the square of the NFW profile, not square the integral, right?
    Am I using the number position distribution and luminosity position distribution in the correct places?
'''

from matplotlib import pyplot as plt
import matplotlib.colors as colors
from math import log
import numpy as np
from matplotlib.lines import Line2D

plt.style.use('latex')


PLOT_SIZE = 50

ALPHA = 1.94
L_MIN_RANGE = [1.0e28, 1.0e34]
L_MAX_RANGE = [1.0e34, 1.0e38]#[1.0e34, 1.0e36]
minPowerStep = (L_MIN_RANGE[1] / L_MIN_RANGE[0]) ** (1.0 / PLOT_SIZE)
maxPowerStep = (L_MAX_RANGE[1] / L_MAX_RANGE[0]) ** (1.0 / PLOT_SIZE)

MULTIPLIER = 5

paperPoint = [1e35, 1e29]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=0.2
TOTAL_FLUX = 7.494712733226778e-10

DRAW_EXTRA_CONTOURS = False
LINE_COLOR = (0.8, 0.3, 0.1)
PATH_TO_FILE = "C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models-position/data-"+ str(MULTIPLIER) + "x/power-law/"
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
    plt.scatter(px, py, marker=('|' if off else '_'), c=LINE_COLOR, sizes = (20,), alpha=0.5)
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
    plt.scatter(px, py, marker=('|' if off else '_'), c=LINE_COLOR, sizes = (20,), alpha=0.5)

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

totalNum = np.stack(totalNum, axis=0)
numSeen = np.stack(numSeen, axis=0)
lumSeen = np.stack(lumSeen, axis=0)

# ========================= Display data =========================

dimMin = len(totalNum)
dimMax = len(totalNum[0])

lMinVals = [L_MIN_RANGE[0] * minPowerStep**i for i in range(PLOT_SIZE)]
lMaxVals = [L_MAX_RANGE[0] * maxPowerStep**j for j in range(PLOT_SIZE)]

fig, ax = plt.subplots(figsize=(6, 4))

plt.xscale("log")
plt.yscale("log")
plt.xlabel("$L_{max}$")
plt.ylabel("$L_{min}$")
plt.title("Power-law, position-dependent{}".format("" if MULTIPLIER is None else (", sensitivity x"+str(MULTIPLIER))))

c1 = plt.pcolor(lMaxVals, lMinVals, totalNum, 
                   norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                   vmax=max([max(v) for v in totalNum])), cmap='Greys_r')
cbar = plt.colorbar(c1, extend='max')
cbar.set_label("$N$")

# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(lMaxVals, lMinVals, numSeen, [10*i for i in range(1, 20)], 
        colors=[(0, i/20.0, 0, 1) for i in range(1, 20)], linewidths=1)
plt.contour(lMaxVals, lMinVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[LINE_COLOR], linewidths=2, label="$N_r=47$")

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(lMaxVals, lMinVals, lumSeen, [0.1*i for i in range(1, 15)], 
        colors=[(1, i/15.0, 1-i/15.0, 1) for i in range(1, 15)], linewidths=1)
plt.contour(lMaxVals, lMinVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed', linewidths=2, label="$R_r=0.2$")

# Observation
shade(numSeen, NUM_PULSARS_ABOVE_THRESHOLD, lMaxVals, lMinVals)
shade(lumSeen, FRAC_ABOVE_THRESHOLD, lMaxVals, lMinVals, True)


# Final points

plt.scatter(paperPoint[0], paperPoint[1], c='purple')

custom_lines = [Line2D([0], [0], color=LINE_COLOR, lw=2),
                Line2D([0], [0], color=LINE_COLOR, lw=2, dashes=(4, 2))]
plt.legend(custom_lines, ["$N_r=47$", "$R_r=0.2$"])
plt.xlim(lMaxVals[0], lMaxVals[-1])
plt.ylim(lMinVals[0], lMinVals[-1])
plt.tight_layout()


if(DRAW_EXTRA_CONTOURS):
    plt.savefig(PATH_TO_FILE + "overlay-extra.png")
if(not DRAW_EXTRA_CONTOURS):
    plt.savefig(PATH_TO_FILE + "power-law-pos-x" + str(MULTIPLIER) + ".png")
    plt.savefig("C:/Users/goods/Dropbox (MIT)/GCE UROP/summaries/jan-2021/figs/power-law/power-law-pos-x" + str(MULTIPLIER) + ".png")

plt.show()