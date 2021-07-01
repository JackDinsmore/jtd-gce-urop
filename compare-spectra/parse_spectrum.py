import numpy as np
from scipy.optimize import curve_fit
import lmfit, os

import warnings
warnings.filterwarnings("ignore")

FERMILAB_GNFW = 0
CALORE = 1
DI_MAURO = 2
AJELLO = 3
FERMILAB_NFW = 4
ABAZAJIAN = 5
GORDON = 6

MY_ROI_SIZE = 0.428821318754
BIG_ROI_SIZE = 0.477550208734
LM_FIT_SCALE = 1e9
ERGS_PER_GEV = 0.00160218
LARGE_ROI_FACTOR = 1.8248083080193496 # Compare 40x40 degree region to 40x40 with mask
GORDON_ROI_FACTOR = 0.7432172914343479 # Compare 7x7 with no mask to 40x40 with mask, with gamma=1.2
ABAZAJIAN_ROI_FACTOR = 0.2912979297594549 # Compare 7x7 with no mask to 40x40 with mask, with gamma=1.1
AJELLO_ROI_FACTOR = 1.2 # Compare 15x15 with no mask to 40x40 with mask, with gamma=1.2

CALORE_L_BREAK, CALORE_ALPHA_BELOW, CALORE_ALPHA_ABOVE = 2.06 * ERGS_PER_GEV, 1.42, 2.63

def calore_power_law(x, scale):
    alpha = np.full_like(x, CALORE_ALPHA_ABOVE)
    alpha[(10**x) < CALORE_L_BREAK] = CALORE_ALPHA_BELOW
    return (10**x)**2 * scale * ((10**x) / CALORE_L_BREAK)**(-alpha)

def power_law(x, scale, alpha_below, alpha_above, l_break):
    alpha = np.full_like(x, alpha_above)
    alpha[10**x < l_break] = alpha_below
    return (10**x)**2 * scale * (10**x / l_break)**(-alpha) * LM_FIT_SCALE

_files = {
    FERMILAB_GNFW:"fermilab-gnfw",
    CALORE:"calore",
    DI_MAURO:"di-mauro",
    AJELLO:"ajello",
    FERMILAB_NFW:"fermilab-nfw",
    ABAZAJIAN:"abazajian",
    GORDON:"gordon"
}
_quadrature_unc = {
    FERMILAB_GNFW: False,
    CALORE: True,
    DI_MAURO: False,
    AJELLO: False,
    FERMILAB_NFW: False,
    ABAZAJIAN: False,
    GORDON: True,
}

class Spectrum:
    def __init__(self, id):
        self.id = id
        self.calore_norm = None
        self.fit_norm = None
        self.fit_alpha_below = None
        self.fit_alpha_above = None
        self.fit_l_break = None

        # Store x in log space, y in linear space.
        # Both CSVs are assumed to be stored in log space
        self.pointsx, self.pointsy = self._load_csv("spectrum-data/"+_files[id]+'/points.csv')
        if os.path.exists("spectrum-data/"+_files[id]+'/up-bands.csv')\
            and os.path.exists("spectrum-data/"+_files[id]+'/down-bands.csv'):
            up_bandsx, up_bandsy = self._load_csv("spectrum-data/"+_files[id]+'/up-bands.csv')
            down_bandsx, down_bandsy = self._load_csv("spectrum-data/"+_files[id]+'/down-bands.csv')
            if _quadrature_unc[id]:
                up_barsx, up_barsy = self._load_csv("spectrum-data/"+_files[id]+'/up-bars.csv')
                down_barsx, down_barsy = self._load_csv("spectrum-data/"+_files[id]+'/down-bars.csv')
                xmin = max(np.min(up_barsx), np.min(up_bandsx), np.min(down_barsx), np.min(down_bandsx))
                xmax = min(np.max(up_barsx), np.max(up_bandsx), np.max(down_barsx), np.max(down_bandsx))
                self.up_barsx = []
                self.up_barsy = []
                self.down_barsx = []
                self.down_barsy = []
                for i, pointx in enumerate(self.pointsx):
                    if pointx < xmin: continue
                    if pointx > xmax: break
                    up_bar = self._get_point(pointx, up_barsx, up_barsy)
                    up_band = self._get_point(pointx, up_bandsx, up_bandsy)
                    down_bar = self._get_point(pointx, down_barsx, down_barsy)
                    down_band = self._get_point(pointx, down_bandsx, down_bandsy)
                    self.up_barsx.append(pointx)
                    self.down_barsx.append(pointx)
                    self.up_barsy.append(np.sqrt((up_bar - self.pointsy[i])**2 + (up_band - self.pointsy[i])**2) + self.pointsy[i])
                    self.down_barsy.append(-np.sqrt((down_bar - self.pointsy[i])**2 + (down_band - self.pointsy[i])**2) + self.pointsy[i])

                self.up_barsx = np.asarray(self.up_barsx)
                self.up_barsy = np.asarray(self.up_barsy)
                self.down_barsx = np.asarray(self.down_barsx)
                self.down_barsy = np.asarray(self.down_barsy)
                print("Used bands and bars for {}, added in quadrature".format(self.get_name()))
            else:
                self.up_barsx = up_bandsx
                self.up_barsy = up_bandsy
                self.down_barsx = down_bandsx
                self.down_barsy = down_bandsy
                print("Used bands for {}".format(self.get_name()))
        else: # There are no bands
            self.up_barsx, self.up_barsy = self._load_csv("spectrum-data/"+_files[id]+'/up-bars.csv')
            self.down_barsx, self.down_barsy = self._load_csv("spectrum-data/"+_files[id]+'/down-bars.csv')
            print("Used bars for {}".format(self.get_name()))

        if id in [DI_MAURO, AJELLO]: # Convert MeV to GeV for certain plots
            self.pointsy /= 1000
            self.up_barsy /= 1000
            self.down_barsy /= 1000

        if id in [FERMILAB_NFW, FERMILAB_GNFW]: # Convert from 2 sigma to 1 sigma
            bars_to_delete = []
            self.up_barsx = list(self.up_barsx)
            self.up_barsy = list(self.up_barsy)
            self.down_barsx = list(self.down_barsx)
            self.down_barsy = list(self.down_barsy)
            for i, barx in enumerate(self.up_barsx):
                point = self._get_point(barx, self.pointsx, self.pointsy)
                if point is None:
                    bars_to_delete.insert(0, i)
                    continue
                self.up_barsy[i] = (self.up_barsy[i] - point) / 2 + point
            # bars_to_delete is sorted from high to low, so this is ok
            for i in bars_to_delete:
                del self.up_barsx[i]
                del self.up_barsy[i]
            bars_to_delete = []
            for i, barx in enumerate(self.down_barsx):
                point = self._get_point(barx, self.pointsx, self.pointsy)
                if point is None:
                    bars_to_delete.insert(0, i)
                    continue
                self.down_barsy[i] = (self.down_barsy[i] - point) / 2 + point
            for i in bars_to_delete:
                del self.down_barsx[i]
                del self.down_barsy[i]
            self.up_barsx = np.asarray(self.up_barsx)
            self.up_barsy = np.asarray(self.up_barsy)
            self.down_barsx = np.asarray(self.down_barsx)
            self.down_barsy = np.asarray(self.down_barsy)

        if id in [ABAZAJIAN]:
            self.pointsy *= 1 / ABAZAJIAN_ROI_FACTOR
            self.up_barsy *= 1 / ABAZAJIAN_ROI_FACTOR
            self.down_barsy *= 1 / ABAZAJIAN_ROI_FACTOR

        if id in [GORDON]:
            self.pointsy *= 1 / GORDON_ROI_FACTOR
            self.up_barsy *= 1 / GORDON_ROI_FACTOR
            self.down_barsy *= 1 / GORDON_ROI_FACTOR

        if id in [AJELLO]:
            self.pointsy *= 1 / AJELLO_ROI_FACTOR
            self.up_barsy *= 1 / AJELLO_ROI_FACTOR
            self.down_barsy *= 1 / AJELLO_ROI_FACTOR

        if id in [DI_MAURO]: # Convert from square ROI to ROI without band
            self.pointsy /= LARGE_ROI_FACTOR
            self.up_barsy /= LARGE_ROI_FACTOR
            self.down_barsy /= LARGE_ROI_FACTOR


    def get_name(self):
        if self.id == FERMILAB_GNFW:
            return "Ref. [14], $\gamma$=1.2"
        elif self.id == CALORE:
            return "Ref. [6]"
        elif self.id == DI_MAURO:
            return "Ref. [7]"
        elif self.id == AJELLO:
            return "Ref. [4]"
        elif self.id == FERMILAB_NFW:
            return "Ref. [14], $\gamma$=1.0"
        elif self.id == ABAZAJIAN:
            return "Ref. [1]"
        elif self.id == GORDON:
            return "Ref. [10]"
        return ""

    def _load_csv(self, filepath):
        # Return x and y values of the spectrum. Both are in logarithmic units. Y is E^2 dN/dE.
        xres = []
        yres = []
        f = open(filepath, 'r')
        for line in f.readlines():
            if line == '': continue
            x, y = line.split(",")
            xres.append(float(x) + np.log10(ERGS_PER_GEV))
            yres.append(10**float(y) * ERGS_PER_GEV)
        f.close()

        # Return the result, scaled to the ROI in question
        if self.id in [ABAZAJIAN, AJELLO, GORDON]:
            return np.asarray(xres), np.asarray(yres)
        if self.id in [DI_MAURO]:
            return np.asarray(xres), np.asarray(yres) * BIG_ROI_SIZE
        if self.id in [FERMILAB_NFW, FERMILAB_GNFW, CALORE]:
            return np.asarray(xres), np.asarray(yres) * MY_ROI_SIZE
        return None

    def _get_point(self, logx, datax, datay):
        i = 0
        if logx < datax[0] or logx > datax[-1]:
            #print("{} is out of range of this function (range {} to {}).".format(logx, datax[0], datax[-1]))
            return
        while datax[i+1] < logx:
            i += 1
        factor = (logx - datax[i]) / (datax[i+1] - datax[i])
        return factor * (datay[i+1] - datay[i]) + datay[i]

    def get_point_log(self, logx):
        return self._get_point(logx, self.pointsx, self.pointsy)

    def get_up_bar_log(self, logx):
        return self._get_point(logx, self.up_barsx, self.up_barsy)

    def get_down_bar_log(self, logx):
        return self._get_point(logx, self.down_barsx, self.down_barsy)

    def _fit_power_law(self, l_min, l_max):
        if self.fit_norm != None:
            return
        log_min = np.log10(l_min) + np.log10(ERGS_PER_GEV)
        log_max = np.log10(l_max) + np.log10(ERGS_PER_GEV)
        model = lmfit.Model(power_law)
        params = lmfit.Parameters()
        params.add(lmfit.Parameter("scale", value=1))
        params.add(lmfit.Parameter("alpha_below", value=CALORE_ALPHA_BELOW))
        params.add(lmfit.Parameter("alpha_above", value=CALORE_ALPHA_ABOVE))
        params.add(lmfit.Parameter("l_break", value=CALORE_L_BREAK))
        this_pointsx = []
        this_pointsy = []
        for i in range(len(self.pointsx)):
            if log_min < self.pointsx[i] < log_max:
                if self.pointsx[i] < np.min(self.up_barsx) or self.pointsx[i] < np.min(self.down_barsx):
                    continue
                if self.pointsx[i] > np.max(self.up_barsx) or self.pointsx[i] > np.max(self.down_barsx):
                    continue
                this_pointsx.append(self.pointsx[i])
                this_pointsy.append(self.pointsy[i] * LM_FIT_SCALE)

        this_pointsx = np.asarray(this_pointsx)
        this_pointsy = np.asarray(this_pointsy)
        bar_widths = [0.5 * (self.get_up_bar_log(l) + self.get_down_bar_log(l)) * LM_FIT_SCALE for l in this_pointsx]
        weights = np.asarray([1.0 / w for w in bar_widths])
        try:
            result = model.fit(data=this_pointsy, params=params, x=this_pointsx, weights=weights)
            self.fit_norm = result.params["scale"].value
            self.fit_alpha_below = result.params["alpha_below"].value
            self.fit_alpha_above = result.params["alpha_above"].value
            self.fit_l_break = result.params["l_break"].value
        except:
            self.fit_norm = -1

    def _fit_calore(self, l_min, l_max):
        if self.calore_norm != None:
            return
        popt, pcov = curve_fit(calore_power_law, self.pointsx, self.pointsy*1000)

        self.calore_norm = popt[0] / 1000

    def get_calore_flux(self, l_min, l_max):
        self._fit_calore(l_min, l_max)
        ergs_l_min = l_min * ERGS_PER_GEV
        ergs_l_max = l_max * ERGS_PER_GEV
        return self.calore_norm * (CALORE_L_BREAK**CALORE_ALPHA_BELOW / (CALORE_ALPHA_BELOW - 2) * (ergs_l_min**(2 - CALORE_ALPHA_BELOW) - CALORE_L_BREAK**(2 - CALORE_ALPHA_BELOW)) +
                            CALORE_L_BREAK**CALORE_ALPHA_ABOVE / (CALORE_ALPHA_ABOVE - 2) * (CALORE_L_BREAK**(2 - CALORE_ALPHA_ABOVE) - ergs_l_max**(2 - CALORE_ALPHA_ABOVE)))

    def get_power_law_flux(self, l_min, l_max):
        self._fit_power_law(l_min, l_max)
        if self.fit_norm < 0:
            return None
        ergs_l_min = l_min * ERGS_PER_GEV
        ergs_l_max = l_max * ERGS_PER_GEV
        return self.fit_norm * (self.fit_l_break**self.fit_alpha_below / (self.fit_alpha_below - 2) * (ergs_l_min**(2 - self.fit_alpha_below) - self.fit_l_break**(2 - self.fit_alpha_below)) +
                            self.fit_l_break**self.fit_alpha_above / (self.fit_alpha_above - 2) * (self.fit_l_break**(2 - self.fit_alpha_above) - ergs_l_max**(2 - self.fit_alpha_above)))

    def get_numerical_flux(self, l_min, l_max):
        log_min = np.log10(l_min * ERGS_PER_GEV)
        log_max = np.log10(l_max * ERGS_PER_GEV)
        if log_min < self.pointsx[0]:
            # Treat all points outside the data as zero flux.
            log_min = self.pointsx[0]
        if log_max > self.pointsx[-1]:
            log_max = self.pointsx[-1]
        # Trapezoidal integration
        integral = 0
        i = 0
        while self.pointsx[i] < log_min:
            i += 1

        integral += 0.5 * (self.pointsy[i] + self.get_point_log(log_min)) * (self.pointsx[i] - log_min)
        while self.pointsx[i+1] < log_max:
            integral += 0.5 * (self.pointsy[i+1] + self.pointsy[i]) * (self.pointsx[i+1] - self.pointsx[i])
            i += 1

        integral += 0.5 * (self.get_point_log(log_max) + self.pointsy[i]) * (log_max - self.pointsx[i])

        return integral * np.log(10)

    def y_data(self, ax, color='k'):
        ax.plot(self.pointsx, self.pointsy, color=color, label=self.get_name())
        thisx = []
        thisy = []
        up_err = []
        down_err = []
        for x in self.pointsx:
            y = self.get_point_log(x)
            u = self.get_up_bar_log(x)
            d = self.get_down_bar_log(x)
            if u is not None and d is not None:
                thisx.append(x)
                thisy.append(y)
                up_err.append(u-y)
                down_err.append(y-d)
        ax.errorbar(thisx, thisy, yerr=[up_err, down_err], color=color, linestyle='none', elinewidth=1)
        ax.set_title(self.get_name())
        ax.set_yscale("log")

    def label_axes(self, ax):
        ax.set_xlabel(self.get_x_label())
        ax.set_ylabel(self.get_y_label())

    def display_data(self, ax, color='k', shape='o'):
        ax.scatter(self.pointsx, self.pointsy, color=color,
            label=self.get_name(), marker=shape)
        thisx = []
        thisy = []
        up_err = []
        down_err = []
        for x in self.pointsx:
            y = self.get_point_log(x)
            u = self.get_up_bar_log(x)
            d = self.get_down_bar_log(x)
            if u is not None and d is not None:
                thisx.append(x)
                thisy.append(y)
                up_err.append(u-y)
                down_err.append(y-d)
        ax.errorbar(thisx, thisy, yerr=[up_err, down_err], color=color,
            linestyle='none', elinewidth=1)
        ax.set_title(self.get_name())
        ax.set_yscale("log")

    def display_calore(self, ax, l_min, l_max, show_all=False):
        self._fit_calore(np.log10(l_min), np.log10(l_max))

        if not show_all:
            calorex = np.linspace(np.log10(l_min * ERGS_PER_GEV),
                np.log10(l_max * ERGS_PER_GEV), 100)
        else:
            calorex = np.linspace(np.min(self.pointsx),
                np.max(self.pointsx), 100)
        calorey = calore_power_law(calorex, self.calore_norm)

        ax.plot(calorex, calorey, color='C1', linestyle='--')

    def display_power_law(self, ax, l_min, l_max, show_all=False):
        self._fit_power_law(l_min, l_max)

        if not show_all:
            fitx = np.linspace(np.log10(l_min * ERGS_PER_GEV),
                np.log10(l_max * ERGS_PER_GEV), 100)
        else:
            fitx = np.linspace(np.min(self.pointsx),
                np.max(self.pointsx), 100)
        fity = power_law(fitx, self.fit_norm, self.fit_alpha_below,
            self.fit_alpha_above, self.fit_l_break) / LM_FIT_SCALE

        ax.plot(fitx, fity, color='C2', linestyle='-.')

    def get_x_label(self):
        return "$\log \\frac{E_\\gamma}{1\\ \\mathrm{erg}}$"

    def get_y_label(self):
        return "$F_\\gamma$ [erg / cm$^2$ / s]"
