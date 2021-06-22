import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model, Parameter, report_fit
import scipy.optimize as optimization

plt.style.use("jcap")

LUM_TO_FLUX = 1.1093417307914119e-46
NUM_POINTS = 100

def power_law(x, c, alpha, lmax):
    return np.log10(c * (10**x)**(-alpha) * np.exp(-10.0**x / lmax))

def log_normal(x, c, l0, sigma):# log(dN/dE)
    return np.log10(c/10**x * np.exp(-(x - np.log10(l0))**2 / (2 * sigma**2)))

def sci_not(i):
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

class Data:
    def __init__(self, name, filename, bar_bottom, bar_top, convert_x_to_log=False,
        convert_y_to_log=False, convert_to_flux=False, convert_to_ergs=False):

        self.name = name
        xs, ys = self.get_data(filename, convert_x_to_log,
            convert_y_to_log, convert_to_flux, convert_to_ergs)
        xdown, ydown = self.get_data(bar_bottom, convert_x_to_log,
            convert_y_to_log, convert_to_flux, convert_to_ergs)
        xup, yup = self.get_data(bar_top, convert_x_to_log,
            convert_y_to_log, convert_to_flux, convert_to_ergs)
        min_x = max(np.min(xs), np.min(xdown), np.min(xup))
        max_x = min(np.max(xs), np.max(xdown), np.max(xup))
        if min_x >= max_x:
            raise Exception("Bars and points share no common domain")
        self.xs = np.linspace(min_x, max_x, 100)
        self.ys = np.asarray([self.at(x, xs, ys) for x in self.xs])
        self.ydown = self.ys - np.asarray([self.at(x, xdown, ydown)
            for x in self.xs])
        self.yup = np.asarray([self.at(x, xup, yup) for x in self.xs]) - self.ys
        self.log_normal_results = None
        self.power_law_results = None

    def plot_data(self, ax):
        ax.scatter(self.xs, self.ys, marker='.')
        ax.errorbar(self.xs, self.ys, yerr=[self.ydown, self.yup])

    def plot_fit(self, ax):
        if self.power_law_results is not None:
            ax.plot(self.xs, log_normal(self.xs, self.power_law_results[0],
            self.power_law_results[2],self.power_law_results[4]),
            label="Power law fit")
        if self.log_normal_results is not None:
            ax.plot(self.xs, log_normal(self.xs, self.log_normal_results[0],
            self.log_normal_results[2],self.log_normal_results[4]),
            label="Log normal fit\n$L_0$={}$\pm${}, $\\sigma$={}$\pm${}".format(
                sci_not(self.log_normal_results[2]),
                sci_not(self.log_normal_results[3]),
                sci_not(self.log_normal_results[4]),
                sci_not(self.log_normal_results[5])))

    def label_axes(self, ax):
        ax.set_xlabel("$\\log(L / $ erg s$^{-1})$")
        ax.set_ylabel("$\\propto L dN / dL$")
        ax.set_title(self.name)

    def plot(self):
        fig, ax = plt.subplots()
        self.plot_data(ax)
        self.plot_fit(ax)
        self.label_axes(ax)
        ax.legend()
        fig.tight_layout()
        return fig

    def fit_power_law(self, c=1, alpha=1, lmax=1e30, display=True):
        model = Model(power_law, independent_vars=['x'])
        result = model.fit(self.ys, x=self.xs,
            c=Parameter('c', value=c),
            alpha=Parameter('alpha', value=alpha),
            lmax=Parameter('lmax', value=lmax),
            weights=1/((self.yup + self.ydown)/2))

        if display:
            report_fit(result.params)

        self.power_law_results = (
            result.params["c"].value, result.params["c"].stderr,
            result.params["alpha"].value, result.params["alpha"].stderr,
            result.params["lmax"].value, result.params["lmax"].stderr)

        return self.power_law_results


    def fit_log_normal(self, c=1e40, l0=1e32, sigma=1, display=True):
        model = model = Model(log_normal, independent_vars=['x'])
        result = model.fit(self.ys, x=self.xs,
            c=Parameter('c', value=c, min=1),
            l0=Parameter('l0', value=l0, min=1),
            sigma=Parameter('sigma', value=sigma, min=1e-3),
            weights=1/((self.yup + self.ydown)/2))

        if display:
            report_fit(result.params)

        self.log_normal_results = (
            result.params["c"].value, result.params["c"].stderr,
            result.params["l0"].value, result.params["l0"].stderr,
            result.params["sigma"].value, result.params["sigma"].stderr)

        return self.log_normal_results


    def get_data(self, filename, convert_x_to_log, convert_y_to_log,
        convert_to_flux, convert_to_ergs):

        xs = []
        ys = []
        f = open(filename, 'r')
        for line in f.readlines():
            if line == '':
                continue
            x, y = line.split(', ')
            x = float(x)
            y = float(y)
            if convert_x_to_log:
                x = np.log10(x)
            if convert_y_to_log:
                y = np.log10(y)
            if convert_to_flux:
                x -= np.log10(LUM_TO_FLUX)
            if convert_to_ergs:
                raise Exception("Unimplemented")
            xs.append(x)
            ys.append(y)
        return np.asarray(xs), np.asarray(ys)

    def at(self, x, xs, ys):
        if x < np.min(xs):
            return 0
        if x > np.max(xs):
            return 0
        i = 0
        while xs[i+1] < x:
            i += 1

        return (ys[i+1] - ys[i]) * (x - xs[i]) / (xs[i+1] - xs[i]) + ys[i]

d = Data("Gautam", "gautam/orange.csv", "gautam/bottom-band.csv", "gautam/top-band.csv",
    convert_to_flux=True)
print(d.fit_log_normal())
d.plot()
plt.show()
