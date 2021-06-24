from matplotlib import pyplot as plt
from math import log, exp
from scipy.special import gammainc, gamma
import matplotlib.colors as colors
from matplotlib.lines import Line2D
import numpy as np

plt.style.use('jcap')


FERMILAB_GNFW_FLUX = 5.716826056794996e-10
CALORE_FLUX = 1.3592975571782718e-09
DI_MAURO_FLUX = 1.2953417255755896e-09
AJELLO_FLUX = 1.8268259089431682e-09
FERMILAB_NFW_FLUX = 7.770902203635269e-10
ABAZAJIAN_FLUX = 4.887121798902831e-10

PLOT_SIZE = 50

L_MIN = 1e29
ALPHA_RANGE = [1.1, 2.5]
L_MAX_RANGE = [1.0e34, 1.0e38]#[1.0e34, 1.0e36]
maxPowerStep = (L_MAX_RANGE[1] / L_MAX_RANGE[0]) ** (1.0 / PLOT_SIZE)

MULTIPLIER = 1

paperPoint = [1e35, 1e29]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=0.2
TOTAL_FLUX =1.2953417255755896e-09

PATH_TO_FILE = "/home/jtdinsmo/Dropbox (MIT)/GCE UROP/luminosity-models-position/data-1x/power-law-alpha/"
SHADE_SCALE=25
LINE_COLOR = "C2"
STYLES = ["solid", "dashed", "dotted", "dashdot"]
SHOW_NUMBERS = False
DRAW_EXTRA_CONTOURS = True

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
    plt.scatter(px, py, marker=('|' if off else '_'), c=LINE_COLOR, sizes = (20,), alpha=1)

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

alphaVals = [ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * (float(i) / PLOT_SIZE) for i in range(PLOT_SIZE)]
lMaxVals = [L_MAX_RANGE[0] * maxPowerStep**j for j in range(PLOT_SIZE)]

fig, ax = plt.subplots(figsize=(6, 4))

#plt.xscale("log")
plt.yscale("log")
plt.xlabel("$\\alpha$")
plt.ylabel("$L_\\mathrm{max}$ [erg / s]")

if SHOW_NUMBERS:
    c1 = plt.pcolor(alphaVals, lMaxVals, totalNum,
                    norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                    vmax=max([max(v) for v in totalNum])), cmap='Greys_r')
    cbar = plt.colorbar(c1, extend='max')
    cbar.set_label("$N_\\text{GCE}$")

# Greens
if DRAW_EXTRA_CONTOURS:
    plt.contour(alphaVals, lMaxVals, numSeen * FERMILAB_GNFW_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=['k'], linewidths=1, linestyles=STYLES[1])
    plt.contour(alphaVals, lMaxVals, numSeen * ABAZAJIAN_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=['k'], linewidths=1, linestyles=STYLES[2])
    plt.contour(alphaVals, lMaxVals, numSeen * AJELLO_FLUX / DI_MAURO_FLUX, [NUM_PULSARS_ABOVE_THRESHOLD], colors=['k'], linewidths=1, linestyles=STYLES[3])
plt.contour(alphaVals, lMaxVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[LINE_COLOR], linewidths=2, linestyles=STYLES[0])

# Reds
plt.contour(alphaVals, lMaxVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed', linewidths=2)

# Observation
if DRAW_EXTRA_CONTOURS:
    shade(numSeen * AJELLO_FLUX / DI_MAURO_FLUX, NUM_PULSARS_ABOVE_THRESHOLD, alphaVals, lMaxVals)
else:
    shade(numSeen, NUM_PULSARS_ABOVE_THRESHOLD, alphaVals, lMaxVals)
shade(lumSeen, FRAC_ABOVE_THRESHOLD, alphaVals, lMaxVals, True)


# Final points

plt.plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='o', markersize=6)

custom_lines = [Line2D([0], [0], color=LINE_COLOR, linestyle=STYLES[0], lw=2),
                Line2D([0], [0], color='k', linestyle=STYLES[1], lw=1),
                Line2D([0], [0], color='k', linestyle=STYLES[2], lw=1),
                Line2D([0], [0], color='k', linestyle=STYLES[3], lw=1),
                Line2D([0], [0], color=LINE_COLOR, linestyle='dashed', lw=2, dashes=(4, 2)),
                Line2D([], [], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='o', linestyle='none', markersize=6),]
plt.legend(custom_lines, ['$N_r$ Di Mauro', "$N_r$ Fermilab gNFW", "$N_r$ Abazajian", "$N_r$ Ajello", '$R_r$', "Fermilab"], loc="lower right")
plt.ylim(lMaxVals[0], lMaxVals[-1])
plt.xlim(alphaVals[0], alphaVals[-1])
plt.tight_layout()


plt.savefig("overlay-power-law-alpha.pdf")

plt.show()
