import numpy as np
import matplotlib.pyplot as plt
import lmfit
import scipy.optimize as optimization

def get_data(filename):
    xs = []
    ys = []
    f = open(filename, 'r')
    for line in f.readlines():
        if line == '':
            continue
        x, y = line.split(', ')
        xs.append(float(x))
        ys.append(float(y))
    return np.asarray(xs), np.asarray(ys)

def power_law(x, c, alpha, lmax):
    return np.log10(c * (10**x)**(-alpha) * np.exp(-10**x / lmax))

def log_normal(x, c, l0, sigma):
    return np.log10(c/10**x * np.exp(-(x - np.log10(l0))**2 / (2 * sigma**2)))

def fit_power_law(filename, color):
    datax, datay = get_data(filename)
    params = optimization.curve_fit(power_law, datax, datay, (1, 1, 1))[0]
    fity = power_law(datax, params[0], params[1], params[2])
    plt.plot(datax, fity, linestyle='dashed', c=color)
    plt.plot(datax, datay, c=color)
    return params[1], params[2]

def fit_log_normal(filename, color):
    datax, datay = get_data(filename)
    params = optimization.curve_fit(log_normal, datax, datay, (1, 1, 1))[0]
    fity = log_normal(datax, params[0], params[1], params[2])
    plt.plot(datax, fity, linestyle='dotted', c=color)
    plt.plot(datax, datay, c=color)
    return params[1], params[2]

plt.plot([], [], linestyle='dashed', c='black', label='power law fit')
plt.plot([], [], linestyle='dotted', c='black', label="log normal fit")
plt.legend()

fit_power_law("blue.csv", color='C0')
fit_power_law("orange.csv", color='C1')
fit_log_normal("blue.csv", color='C0')
fit_log_normal("orange.csv", color='C1')
plt.xlabel("log10(F/(ergs/s/cm^2))")
plt.ylabel("unknown")

plt.savefig("fit results.png")
plt.show()