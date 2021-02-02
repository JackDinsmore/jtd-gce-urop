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

plt.style.use('latex')


POWER_STEP = 1.1 # 1 is the minimum

ALPHA = 1.94
L_MIN_RANGE = [1.0e28, 1.0e34]
L_MAX_RANGE = [1.0e34, 1.0e36]

paperPoint = [1e35, 1e29]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=0.2

DRAW_EXTRA_CONTOURS = False
PATH_TO_FILE = "C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models-position/data/power-law/"

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

totalNum = np.stack(totalNum, axis=0)
numSeen = np.stack(numSeen, axis=0)
lumSeen = np.stack(lumSeen, axis=0)

# ========================= Display data =========================

dimMin = len(totalNum)
dimMax = len(totalNum[0])

lMinVals = [L_MIN_RANGE[0] * POWER_STEP**i for i in range(dimMin)]
lMaxVals = [L_MAX_RANGE[0] * POWER_STEP**j for j in range(dimMax)]

fig, ax = plt.subplots(figsize=(6,4))
plt.text(0.95, 0.95, 'Green: number limit\nRed: luminosity limit', 
    horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, color='white', backgroundcolor=(0, 0, 0, 0.3))

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Lmax")
plt.ylabel("Lmin")
plt.title("Exp cutoff (alpha={0})".format(ALPHA))

c1 = plt.pcolor(lMaxVals, lMinVals, totalNum, 
                   norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                   vmax=max([max(v) for v in totalNum])), cmap='PuBu_r')
cbar = plt.colorbar(c1, extend='max')
cbar.set_label("Number of MSPs")

# Greens
if(DRAW_EXTRA_CONTOURS):
    plt.contour(lMaxVals, lMinVals, numSeen, [10*i for i in range(1, 20)], 
        colors=[(0, i/20.0, 0, 1) for i in range(1, 20)], linewidths=1)
plt.contour(lMaxVals, lMinVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[(0, 1, 0)], linewidths=2)

# Reds
if(DRAW_EXTRA_CONTOURS):
    plt.contour(lMaxVals, lMinVals, lumSeen, [0.1*i for i in range(1, 15)], 
        colors=[(1, i/15.0, 1-i/15.0, 1) for i in range(1, 15)], linewidths=1)
plt.contour(lMaxVals, lMinVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[(1, 0, 0)], linewidths=2)


plt.scatter(paperPoint[0], paperPoint[1], c='purple')


if(DRAW_EXTRA_CONTOURS):
    plt.savefig(PATH_TO_FILE+"overlay-extra.png")
if(not DRAW_EXTRA_CONTOURS):
    plt.savefig(PATH_TO_FILE+"overlay.png")

plt.show()