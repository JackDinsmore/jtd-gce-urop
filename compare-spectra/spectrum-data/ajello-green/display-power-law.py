import matplotlib.pyplot as plt
import numpy as np

NORM = 1.4e-4
L_MAX = 3.0
ALPHA = 0.4

def power_law(x, norm, lmax, alpha):
    return norm * x**2 * x**(-alpha) * np.exp(-x / lmax)

def load_csv(filepath):
    # Return x and y values of the spectrum. Both are in logarithmic units. Y is E^2 dN/dE.
    xres = []
    yres = []
    f = open(filepath, 'r')
    for line in f.readlines():
        if line == '': continue
        x, y = line.split(",")
        xres.append(10**float(x))
        yres.append(10**float(y))
    f.close()

    return np.asarray(xres), np.asarray(yres)

def get_point(x, datax, datay):
    i = 0
    if x < datax[0] or x > datax[-1]:
        #print("{} is out of range of this function (range {} to {}).".format(logx, datax[0], datax[-1]))
        return
    while datax[i+1] < x:
        i += 1
    factor = (x - datax[i]) / (datax[i+1] - datax[i])
    return factor * (datay[i+1] - datay[i]) + datay[i]

downx, downy = load_csv("down-bands.csv")
upx, upy = load_csv("up-bands.csv")
ptsx = 10**np.linspace(np.log10(np.min(downx)), np.log10(np.max(downx)), 10)
ptsy = power_law(ptsx, NORM, L_MAX, ALPHA)
fake_ptsy = [(get_point(x, downx, downy) + get_point(x, upx, upy)) / 2 for x in ptsx]

plt.plot(downx, downy, color='k')
plt.plot(upx, upy, color='k')
plt.plot(ptsx, ptsy, label="Fit")
plt.plot(ptsx, fake_ptsy, label="Mid points")

f = open("points.csv", 'w')
for i in range(len(ptsx)):
    f.write("{}, {}\n".format(np.log10(ptsx[i]), np.log10(ptsy[i])))
f.close()

plt.xscale("log")
plt.yscale("log")
plt.xlim(1, 100)
plt.ylim(7e-6, 6e-4)
plt.xlabel("E (GeV)")
plt.ylabel("E^2 dN/dE (MeV / cm^2 / s)")
plt.legend()

plt.show()
