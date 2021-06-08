import numpy as np
from scipy.optimize import curve_fit
import lmfit

import warnings
warnings.filterwarnings("ignore")

FERMILAB_GNFW = 0
CALORE = 1
DI_MAURO = 2
AJELLO = 3
FERMILAB_NFW = 4
ABAZAJIAN = 5
GORDON = 6

BANDS_OR_BARS = 'bands'# 'bars'

ROI_SIZE = 0.428821318754
LM_FIT_SCALE = 1e9
ERGS_PER_GEV = 0.00160218
LARGE_ROI_FACTOR = 1.8248083080193496#2.101243483243096 # Compare 40x40 degree region to 40x40 with mask
GORDON_ROI_FACTOR = 0.7432172914343479#1.0138798168363303 # Compare 7x7 with no mask to 40x40 with mask, with gamma=1.2
ABAZAJIAN_ROI_FACTOR = 0.2912979297594549#0.4208269993877883 # Compare 7x7 with no mask to 40x40 with mask, with gamma=1

CALORE_L_BREAK, CALORE_ALPHA_BELOW, CALORE_ALPHA_ABOVE = 2.06 * ERGS_PER_GEV, 1.42, 2.63

def calore_power_law(x, scale):
    alpha = np.full_like(x, CALORE_ALPHA_ABOVE)
    alpha[(10**x) < CALORE_L_BREAK] = CALORE_ALPHA_BELOW
    return (10**x)**2 * scale * ((10**x) / CALORE_L_BREAK)**(-alpha)

def power_law(x, scale, alpha_below, alpha_above, l_break):
    alpha = np.full_like(x, alpha_above)
    alpha[10**x < l_break] = alpha_below
    return (10**x)**2 * scale * (10**x / l_break)**(-alpha) * LM_FIT_SCALE

_files = { FERMILAB_GNFW:"fermilab-gnfw", CALORE:"calore", DI_MAURO:"di-mauro", AJELLO:"ajello", FERMILAB_NFW:"fermilab-nfw", ABAZAJIAN:"abazajian", GORDON:"gordon" }

class Spectrum:
    def __init__(self, id):
        self.pointsx, self.pointsy = self._load_csv("spectrum-data/"+_files[id]+'/points.csv')
        try:
            self.up_barsx, self.up_barsy = self._load_csv("spectrum-data/"+_files[id]+'/up-'+BANDS_OR_BARS + '.csv')
            self.down_barsx, self.down_barsy = self._load_csv("spectrum-data/"+_files[id]+'/down-'+BANDS_OR_BARS + '.csv')
        except: # There are no bands
            self.up_barsx, self.up_barsy = self._load_csv("spectrum-data/"+_files[id]+'/up-bars.csv')
            self.down_barsx, self.down_barsy = self._load_csv("spectrum-data/"+_files[id]+'/down-bars.csv')


        if id in [AJELLO, DI_MAURO]: # Convert MeV to GeV for certain plots
            self.pointsy /= 1000
            self.up_barsy /= 1000
            self.down_barsy /= 1000

        if id in [ABAZAJIAN]:
            self.pointsy *= 1 / ABAZAJIAN_ROI_FACTOR
            self.up_barsy *= 1 / ABAZAJIAN_ROI_FACTOR
            self.down_barsy *= 1 / ABAZAJIAN_ROI_FACTOR

        if id in [GORDON]:
            self.pointsy *= 1 / GORDON_ROI_FACTOR
            self.up_barsy *= 1 / GORDON_ROI_FACTOR
            self.down_barsy *= 1 / GORDON_ROI_FACTOR

        if id in [DI_MAURO, AJELLO]: # Convert from square ROI to ROI without band
            self.pointsy /= LARGE_ROI_FACTOR
            self.up_barsy /= LARGE_ROI_FACTOR
            self.down_barsy /= LARGE_ROI_FACTOR
        
        self.id = id
        self.calore_norm = None
        self.fit_norm = None
        self.fit_alpha_below = None
        self.fit_alpha_above = None
        self.fit_l_break = None
    
    def get_name(self):
        if self.id == FERMILAB_GNFW:
            return "Fermilab gNFW"
        elif self.id == CALORE:
            return "Calore"
        elif self.id == DI_MAURO:
            return "Di Mauro"
        elif self.id == AJELLO:
            return "Ajello"
        elif self.id == FERMILAB_NFW:
            return "Fermilab NFW"
        elif self.id == ABAZAJIAN:
            return "Abazajian"
        elif self.id == GORDON:
            return "Gordon"
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
        return np.asarray(xres), np.asarray(yres) * ROI_SIZE
    
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
        if log_min < self.pointsx[0] or log_max > self.pointsx[-1]:
            print("Integration out of range")
            print(self.pointsx[0], self.pointsx[-1])
            print(log_min, log_max)
            return

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

    def display_data(self, ax, color='k'):
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
        ax.errorbar(thisx, thisy, yerr=[up_err, down_err], color=color, linestyle='none')
        ax.set_title(self.get_name())
        ax.set_yscale("log")

    def label_axes(self, ax):
        ax.set_xlabel(self.get_x_label())
        ax.set_ylabel(self.get_y_label())

    def display_calore(self, ax, l_min, l_max):
        self._fit_calore(np.log10(l_min), np.log10(l_max))

        calorex = np.linspace(np.log10(l_min * ERGS_PER_GEV), np.log10(l_max * ERGS_PER_GEV), 100)
        calorey = calore_power_law(calorex, self.calore_norm)
        
        ax.plot(calorex, calorey, color='C1', linestyle='--')

    def display_power_law(self, ax, l_min, l_max):
        self._fit_power_law(l_min, l_max)

        if self.fit_norm < 1:
            return

        fitx = np.linspace(np.log10(l_min * ERGS_PER_GEV), np.log10(l_max * ERGS_PER_GEV), 100)
        fity = power_law(fitx, self.fit_norm, self.fit_alpha_below, self.fit_alpha_above, self.fit_l_break) / LM_FIT_SCALE
        
        ax.plot(fitx, fity, color='C2', linestyle='--')

    def get_x_label(self):
        return "$\log \\frac{E}{1\\ \\mathrm{erg}}$"

    def get_y_label(self):
        return "$E^2 \\frac{dN}{dE}$ [erg / cm$^2$ / s]"