import numpy as np
import parse_spectrum as ps
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib.transforms as mtrans

plt.style.use("jcap")

NAME = ps.DI_MAURO
SPECTRUM_RANGE = [0.5, 10]# GeV
NUM_PLOTS = 7
SHAPES=['o', 's', '^', '*', 'd', '+', 'x']

fig, ax = plt.subplots()

f = ps.Spectrum(NAME)
f.display_data(ax)
f.display_calore(ax, SPECTRUM_RANGE[0], SPECTRUM_RANGE[1], show_all=True)
f.display_power_law(ax, SPECTRUM_RANGE[0], SPECTRUM_RANGE[1], show_all=True)
numerical_flux = f.get_numerical_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
calore_flux = f.get_calore_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
fit_flux = f.get_power_law_flux(SPECTRUM_RANGE[0], SPECTRUM_RANGE[1])
f.label_axes(ax)

ax.plot([], [], color='C1', linestyle='--', label="Calore fit")
ax.plot([], [], color='C2', linestyle='-.', label="Power law fit")

plt.title("")
fig.tight_layout()
ax.legend()
fig.savefig("example.pdf")
plt.show()
