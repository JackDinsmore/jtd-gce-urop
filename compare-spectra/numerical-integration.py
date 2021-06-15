import numpy as np
import parse_spectrum as ps
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib.transforms as mtrans

plt.style.use("latex")

SPECTRUM_RANGE = [0.5, 10]# GeV
NUM_PLOTS = 7

fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(10, 6), sharex='col', sharey='row')

allfig, allax = plt.subplots()

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
names = []

for i in range(0, NUM_PLOTS):
    ax = axes[i%3][i//3]
    f = ps.Spectrum(i)
    names.append(f.get_name())
    f.display_data(allax, color="C"+str(i))
    f.display_data(ax)
    f.display_calore(ax, SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    f.display_power_law(ax, SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    numerical_flux = f.get_numerical_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    calore_flux = f.get_calore_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    fit_flux = f.get_power_law_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
    num_fluxes.append(numerical_flux)
    calore_fluxes.append(calore_flux)
    fit_fluxes.append(fit_flux)
    if i//3 == 0:
        ax.set_ylabel(f.get_y_label())
    if i % 3 == 1:
        ax.set_xlabel(f.get_x_label())

    f.label_axes(allax)

    print(f.get_name(), numerical_flux)
    
    ax.annotate("Num: {}\nCalore: {}\nFit: {}".format(sciNot(numerical_flux), sciNot(calore_flux), sciNot(fit_flux)), (0.2, 0.1), xycoords='axes fraction')


axes[(NUM_PLOTS-1)%3][(NUM_PLOTS-1)//3].set_xlabel(f.get_x_label())

fig.delaxes(axes[2][2])
fig.delaxes(axes[1][2])

while None in fit_fluxes:
    index = fit_fluxes.index(None)
    del fit_fluxes[index]
    del calore_fluxes[index]
    del num_fluxes[index]
    del names[index]

while None in num_fluxes:
    index = num_fluxes.index(None)
    del fit_fluxes[index]
    del calore_fluxes[index]
    del num_fluxes[index]
    del names[index]

x = np.arange(len(fit_fluxes))
width = 0.25
hist_fig, hist_ax = plt.subplots()
hist_ax.bar(x, num_fluxes, width, label="Numerical")
hist_ax.bar(x + width, calore_fluxes, width, label="Calore")
hist_ax.bar(x + 2 * width, fit_fluxes, width, label="Power law")
hist_ax.legend()
hist_ax.set_ylabel(f.get_y_label())
hist_ax.set_xticks(x + width)
hist_ax.set_xticklabels(names)
hist_fig.savefig("integral-hist.png")



calore_fit = mlines.Line2D([], [], color='C1', linestyle='--', label='Calore fit')
power_law_fit = mlines.Line2D([], [], color='C2', linestyle='--', label='Power law fit')
fig.legend(handles=[calore_fit, power_law_fit], loc='lower right')


allfig.legend()
allax.set_title("")

fig.tight_layout()
allfig.tight_layout()


fig.savefig("summary.png")
allfig.savefig("all-spectra.png")
plt.show()