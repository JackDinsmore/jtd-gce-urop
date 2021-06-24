from matplotlib import pyplot as plt
import matplotlib.colors as colors
from math import log
import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import ImageGrid

plt.style.use('jcap')


PLOT_SIZE = 50

ALPHA = 1.94
L_MIN_RANGE = [1.0e28, 1.0e34]
L_MAX_RANGE = [1.0e34, 1.0e38]#[1.0e34, 1.0e36]
minPowerStep = (L_MIN_RANGE[1] / L_MIN_RANGE[0]) ** (1.0 / PLOT_SIZE)
maxPowerStep = (L_MAX_RANGE[1] / L_MAX_RANGE[0]) ** (1.0 / PLOT_SIZE)

MULTIPLIER = 1

paperPoint = [1e35, 1e29]

NUM_PULSARS_ABOVE_THRESHOLD = 47
FRAC_ABOVE_THRESHOLD=0.2
TOTAL_FLUX = 1.2953417255755896e-09#7.494712733226778e-10

DRAW_EXTRA_CONTOURS = False
LINE_COLOR ="C0"
SHADE_SCALE=25

def get_path(mult):
    return "/home/jtdinsmo/Dropbox (MIT)/GCE UROP/luminosity-models-position/data-"\
        + str(mult) + "x/power-law/"

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

lMinVals = [L_MIN_RANGE[0] * minPowerStep**i for i in range(PLOT_SIZE)]
lMaxVals = [L_MAX_RANGE[0] * maxPowerStep**j for j in range(PLOT_SIZE)]

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
    ax.set_xlim(left=min(lMaxVals), right=max(lMaxVals))
    ax.set_ylim(bottom=min(lMinVals), top=max(lMinVals))
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel("$L_\\mathrm{min}$ [ergs / s]")
axs[-1].set_xlabel("$L_\\mathrm{max}$ [ergs / s]")

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

    totalNum = np.stack(totalNum, axis=0)
    numSeen = np.stack(numSeen, axis=0)
    lumSeen = np.stack(lumSeen, axis=0)

    # ========================= Display data =========================

    c1 = axs[i].contourf(lMaxVals, lMinVals, totalNum,
                       norm=colors.LogNorm(vmin=min([min(v) for v in totalNum]),
                       vmax=max([max(v) for v in totalNum])), cmap='Greys_r')

    axs[i].contour(lMaxVals, lMinVals, numSeen, [NUM_PULSARS_ABOVE_THRESHOLD], colors=[LINE_COLOR])
    axs[i].contour(lMaxVals, lMinVals, lumSeen, [FRAC_ABOVE_THRESHOLD], colors=[LINE_COLOR], linestyles='dashed')

    shade(axs[i], numSeen, NUM_PULSARS_ABOVE_THRESHOLD, lMaxVals, lMinVals)
    shade(axs[i], lumSeen, FRAC_ABOVE_THRESHOLD, lMaxVals, lMinVals, True)

    axs[i].plot(paperPoint[0], paperPoint[1], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='o', markersize=6)
    i += 1

cbar = axs[-1].cax.colorbar(c1)
cbar.set_label("$N_\\textrm{GCE}$")

custom_lines = [Line2D([0], [0], color=LINE_COLOR),
                Line2D([0], [0], color=LINE_COLOR, dashes=(4, 2)),
                Line2D([], [], markeredgecolor='black', markerfacecolor=LINE_COLOR, marker='o', linestyle='None', markersize=6),]
plt.legend(custom_lines, ["$N_r=47$", "$R_r=0.2$", "Fermilab point"], loc='lower right')
plt.xlim(lMaxVals[0], lMaxVals[-1])
plt.ylim(lMinVals[0], lMinVals[-1])
plt.tight_layout()

plt.savefig("power-law-sensitivity.pdf")

plt.show()
