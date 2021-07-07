import numpy as np
import parse_spectrum as ps
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib.transforms as mtrans

plt.style.use("jcap")


SPECTRUM_RANGE = [0.1, 100]# GeV
SPECTRUM_RANGE_LOW = [0.1, 10]# GeV
NUM_PLOTS = 9
SHAPES = ['o', 's', '*', 'd', 'x', 'o', 's', '<', '>']
FILL_STYLE = [None, None, 'none', None, None, 'none', 'none', 'none', 'none']

fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10, 6), sharex='col', sharey='row')

allfig, allax = plt.subplots(figsize=(6.8, 4.5))


allax.set_ylim(1e-8, 6e-7)

def sciNot(i):
    if i == None:
        return "None"
    if i == np.nan:
        return "NaN"
    if i == 0:
        return "0"
    negative = False
    if i < 0:
        negative = True
        i *= -1
    d = int(np.log10(i))
    if np.log10(i) < 0: d -= 1
    if not negative:
        return str(round(i / 10**d, 3)) + "e"+str(d)
    return "-"+str(round(i / 10**d, 3)) + "e"+str(d)

num_fluxes = []
calore_fluxes = []
fit_fluxes = []
num_fluxes_low = []
calore_fluxes_low = []
fit_fluxes_low = []
names = []

calore_flux_uncs = []
fit_flux_uncs = []
calore_flux_low_uncs = []
fit_flux_low_uncs = []

for i in range(0, NUM_PLOTS):
    ax = axes[i%3][i//3]
    f = ps.Spectrum(i)
    names.append(f.get_name())
    f.display_data(allax, color="C"+str(i), shape=SHAPES[i],
        fill_color=FILL_STYLE[i])
    f.display_data(ax)
    f.display_calore(ax, SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    f.display_power_law(ax, SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    numerical_flux = f.get_numerical_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    calore_flux = f.get_calore_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    fit_flux = f.get_power_law_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    num_fluxes.append(numerical_flux)
    calore_fluxes.append(calore_flux[0])
    fit_fluxes.append(fit_flux[0])
    calore_flux_uncs.append(calore_flux[1])
    fit_flux_uncs.append(fit_flux[1])
    if i//3 == 0:
        ax.set_ylabel(f.get_y_label())
    if i % 3 == 1:
        ax.set_xlabel(f.get_x_label())

    f.label_axes(allax)

    print(f.get_name(), numerical_flux)

    ax.annotate("Num: {}\nCalore: {}\nFit: {}".format(sciNot(numerical_flux), sciNot(calore_flux[0]), sciNot(fit_flux[0])), (0.2, 0.1), xycoords='axes fraction', size=8)


    numerical_flux_low = f.get_numerical_flux(SPECTRUM_RANGE_LOW[0], SPECTRUM_RANGE_LOW[1], override=True)
    calore_flux_low = f.get_calore_flux(SPECTRUM_RANGE_LOW[0], SPECTRUM_RANGE_LOW[1], override=True)
    fit_flux_low = f.get_power_law_flux(SPECTRUM_RANGE_LOW[0], SPECTRUM_RANGE_LOW[1], override=True)
    num_fluxes_low.append(numerical_flux_low)
    calore_fluxes_low.append(calore_flux_low[0])
    fit_fluxes_low.append(fit_flux_low[0])
    calore_flux_low_uncs.append(calore_flux[1])
    fit_flux_low_uncs.append(fit_flux[1])


axes[(NUM_PLOTS-1)%3][(NUM_PLOTS-1)//3].set_xlabel(f.get_x_label())

#fig.delaxes(axes[2][2])

while None in fit_fluxes:
    index = fit_fluxes.index(None)
    fit_fluxes[index] = 0
    #del fit_fluxes[index]
    #del calore_fluxes[index]
    #del num_fluxes[index]
    #del names[index]

x = np.arange(len(fit_fluxes))
width = 0.25
delta = 0#width/6

hist_fig, hist_ax = plt.subplots(figsize=(8, 5))
hist_ax.bar(x, num_fluxes, width, label="Numerical ({}-{} GeV)".format(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1]), color='purple', alpha=0.3)
hist_ax.bar(x + width, calore_fluxes, width, label="Ref [6] fit", linewidth=0.5, alpha=0.5, color='r')
hist_ax.bar(x + 2 * width, fit_fluxes, width, label="Power law fit", linewidth=0.5, alpha=0.5, color='b')
hist_ax.errorbar(x + width + delta, calore_fluxes, yerr=calore_flux_uncs, linestyle="none", linewidth=1, color='r')
hist_ax.errorbar(x + 2 * width + delta, fit_fluxes, yerr=fit_flux_uncs, linestyle="none", linewidth=1, color='b')

hist_ax.bar(x, num_fluxes_low, width, label="{}-{} GeV".format(SPECTRUM_RANGE_LOW[0], SPECTRUM_RANGE_LOW[1]), fill=False, edgecolor='k')
hist_ax.bar(x + width, calore_fluxes_low, width, linewidth=1, fill=False, edgecolor='k')
hist_ax.bar(x + 2 * width, fit_fluxes_low, width, linewidth=1, fill=False, edgecolor='k')
#hist_ax.errorbar(x + width - delta, calore_fluxes_low, yerr=calore_flux_low_uncs, linestyle="none", linewidth=1, color='k')
#hist_ax.errorbar(x + 2 * width - delta, fit_fluxes_low, yerr=fit_flux_low_uncs, linestyle="none", linewidth=1, color='k')

hist_ax.legend()
hist_ax.set_ylabel("$F_\\mathrm{GCE}$ [erg / cm$^2$ / s]")
hist_ax.set_xticks(x + width)
hist_ax.set_xticklabels(names, rotation=70, size=11)
hist_fig.tight_layout()
hist_fig.savefig("integral-hist.pdf")



calore_fit = mlines.Line2D([], [], color='C1', linestyle='--', label='Calore fit')
power_law_fit = mlines.Line2D([], [], color='C2', linestyle='--', label='Power law fit')
fig.legend(handles=[calore_fit, power_law_fit], loc='lower right')


allax.legend(loc='lower left', ncol=3)
allax.set_title("")

fig.tight_layout()
allfig.tight_layout()


fig.savefig("summary.pdf")
allfig.savefig("all-spectra.pdf")
plt.show()
