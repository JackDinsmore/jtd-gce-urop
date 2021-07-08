import numpy as np
from scipy.optimize import curve_fit
import lmfit, os

import warnings
warnings.filterwarnings("ignore")

DI_MAURO = 0
CALORE = 1
FERMILAB_GNFW = 2
FERMILAB_NFW = 3
ABAZAJIAN = 4
GORDON = 5
AJELLO_2017 = 8
AJELLO_BLACK = 6
AJELLO_GREEN = 7
ARTIFICIAL = 9

MY_ROI_SIZE = 0.428821318754
BIG_ROI_SIZE = 0.477550208734
LM_FIT_SCALE = 1e9
ERGS_PER_GEV = 0.00160218
LARGE_ROI_FACTOR = 1.8961920419764398 # Compare 40x40 degree region to 40x40 with mask
GORDON_ROI_FACTOR = 0.8685899550623982 # Compare 7x7 with no mask to 40x40 with mask, with gamma=1.2
ABAZAJIAN_ROI_FACTOR = 0.5496547377000791# Compare 7x7 with no mask to 40x40 with mask, with gamma=1.1 # UNKNOWN
AJELLO_ROI_FACTOR = 1.2624815974139625 # Compare 15x15 with no mask to 40x40 with mask, with gamma=1.2

Y_LIM = 1e-8 * ERGS_PER_GEV# lOWERMOST Y LIMIT
ERROR_BAR_SIZE = 0.2

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
    AJELLO_BLACK:"ajello-black",
    AJELLO_GREEN:"ajello-green",
    AJELLO_2017:"ajello-2017",
    FERMILAB_NFW:"fermilab-nfw",
    ABAZAJIAN:"abazajian",
    GORDON:"gordon",
    ARTIFICIAL:"artificial"
}
_quadrature_unc = {
    FERMILAB_GNFW: False,
    CALORE: True,
    DI_MAURO: False,
    AJELLO_BLACK: False,
    AJELLO_GREEN: False,
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
                xmin_down = max(np.min(down_barsx), np.min(down_bandsx))
                xmax_down = min(np.max(down_barsx), np.max(down_bandsx))
                xmin_up = max(np.min(up_barsx), np.min(up_bandsx))
                xmax_up = min(np.max(up_barsx), np.max(up_bandsx))
                self.up_barsx = []
                self.up_barsy = []
                self.down_barsx = []
                self.down_barsy = []
                for i, pointx in enumerate(self.pointsx):
                    if xmin_down < pointx < xmax_down:
                        down_bar = self._get_point(pointx, down_barsx, down_barsy)
                        down_band = self._get_point(pointx, down_bandsx, down_bandsy)
                        self.down_barsx.append(pointx)
                        self.down_barsy.append(-np.sqrt((down_bar - self.pointsy[i])**2 + (down_band - self.pointsy[i])**2) + self.pointsy[i])
                    if xmin_up < pointx < xmax_up:
                        up_bar = self._get_point(pointx, up_barsx, up_barsy)
                        up_band = self._get_point(pointx, up_bandsx, up_bandsy)
                        self.up_barsx.append(pointx)
                        self.up_barsy.append(np.sqrt((up_bar - self.pointsy[i])**2 + (up_band - self.pointsy[i])**2) + self.pointsy[i])

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

        if id in [DI_MAURO, AJELLO_BLACK, AJELLO_GREEN]: # Convert MeV to GeV for certain plots
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

        if id in [AJELLO_BLACK, AJELLO_GREEN]:
            self.pointsy *= 1 / AJELLO_ROI_FACTOR
            self.up_barsy *= 1 / AJELLO_ROI_FACTOR
            self.down_barsy *= 1 / AJELLO_ROI_FACTOR

        if id in [DI_MAURO]: # Convert from square ROI to ROI without band
            self.pointsy /= LARGE_ROI_FACTOR
            self.up_barsy /= LARGE_ROI_FACTOR
            self.down_barsy /= LARGE_ROI_FACTOR


    def get_name(self):
        if self.id == ARTIFICIAL:
            return "Fake"
        elif self.id == FERMILAB_GNFW:
            return "Zhong 2020, $\gamma$=1.2"
        elif self.id == CALORE:
            return "Calore 2015"
        elif self.id == DI_MAURO:
            return "Di Mauro 2021"
        elif self.id == AJELLO_GREEN:
            return "Ajello 2016, OB"
        elif self.id == AJELLO_BLACK:
            return "Ajello 2016, PSR"
        elif self.id == FERMILAB_NFW:
            return "Zhong 2020, $\gamma$=1.0"
        elif self.id == ABAZAJIAN:
            return "Abazajian 2014"
        elif self.id == GORDON:
            return "Gordon 2013"
        elif self.id == AJELLO_2017:
            return "Ajello 2017"
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
        if self.id in [ABAZAJIAN, AJELLO_BLACK, AJELLO_GREEN, GORDON, ARTIFICIAL, AJELLO_2017]:
            return np.asarray(xres), np.asarray(yres)
        if self.id in [DI_MAURO]:
            return np.asarray(xres), np.asarray(yres) * BIG_ROI_SIZE
        if self.id in [FERMILAB_NFW, FERMILAB_GNFW, CALORE]:
            return np.asarray(xres), np.asarray(yres) * MY_ROI_SIZE
        return None

    def _get_point(self, logx, datax, datay):
        i = 0
        if logx < datax[0] or logx > datax[-1]:
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

    def _fit_power_law(self, l_min, l_max, override):
        if self.fit_norm != None and not override:
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
            self.fit_norm_unc = result.params["scale"].stderr
            self.fit_alpha_below = result.params["alpha_below"].value
            self.fit_alpha_below_unc = result.params["alpha_below"].stderr
            self.fit_alpha_above = result.params["alpha_above"].value
            self.fit_alpha_above_unc = result.params["alpha_above"].stderr
            self.fit_l_break = result.params["l_break"].value
            self.fit_l_break_unc = result.params["l_break"].stderr
        except:
            self.fit_norm = -1

    def _fit_calore(self, l_min, l_max, override):
        """if self.calore_norm != None and not override:
            return
        popt, pcov = curve_fit(calore_power_law, self.pointsx, self.pointsy*1000)

        self.calore_norm = popt[0] / 1000"""
        if self.calore_norm != None and not override:
            return
        log_min = np.log10(l_min) + np.log10(ERGS_PER_GEV)
        log_max = np.log10(l_max) + np.log10(ERGS_PER_GEV)
        model = lmfit.Model(calore_power_law)
        params = lmfit.Parameters()
        params.add(lmfit.Parameter("scale", value=1))
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
            self.calore_norm = result.params["scale"].value / LM_FIT_SCALE
            self.calore_norm_unc = result.params["scale"].stderr / LM_FIT_SCALE
        except:
            self.calore_norm = -1

    def get_calore_flux(self, l_min, l_max, override=False):
        self._fit_calore(l_min, l_max, override)
        a = l_min * ERGS_PER_GEV
        b = l_max * ERGS_PER_GEV
        c = self.calore_norm
        L = CALORE_L_BREAK
        a1 = CALORE_ALPHA_BELOW
        a2 = CALORE_ALPHA_ABOVE
        val = c * ((L**2 - a**(2 - a1) * L**a1) / (2 - a1) - (L**2 - b**(2 - a2) * L**a2) / (2 - a2))
        return (val, val * self.calore_norm_unc / c)

    def get_power_law_flux(self, l_min, l_max, override=False):
        self._fit_power_law(l_min, l_max, override)
        if self.fit_norm < 0:
            return None

        c = self.fit_norm
        L = self.fit_l_break
        a1 = self.fit_alpha_below
        a2 = self.fit_alpha_above
        a = l_min * ERGS_PER_GEV
        b = l_max * ERGS_PER_GEV
        val = c * ((L**2 - a**(2 - a1) * L**a1) / (2 - a1) - (L**2 - b**(2 - a2) * L**a2) / (2 - a2))

        norm_unc = val * self.fit_norm_unc / c
        break_unc = c * self.fit_l_break_unc * ((2 * L - a1 * a**(2 - a1) * L**(a1-1)) / (2 - a1) - (2 * L - b**(2 - a2) * a2 * L**(a2-1)) / (2 - a2))
        below_unc = c * self.fit_alpha_below_unc * a**(-a1) / (a1 - 2)**2 * \
             (a**(a1) * L**2 - a**2 * L**a1 - (a1 - 2)* a**2 * L**a1 * np.log(a) + (a1 -2)*a**2 * L**a1 * np.log(L))
        above_unc = c * self.fit_alpha_above_unc * b**(-a2) / (a2 - 2)**2 * \
            (-b**(a2) * L**2 + b**2 * L**a2 + (a2 - 2)* b**2 * L**a2 * np.log(b) - (a2 -2)*b**2 * L**a2 * np.log(L))
        
        return (val, np.sqrt(norm_unc**2 + break_unc**2 + above_unc**2 + below_unc**2))

    def get_numerical_flux(self, l_min, l_max, override=False):
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

    def label_axes(self, ax):
        ax.set_xlabel(self.get_x_label())
        ax.set_ylabel(self.get_y_label())

    def display_data(self, ax, color='k', shape='o', fill_color=None):
        fcolor = color if None else fill_color
        ax.scatter(10**self.pointsx / ERGS_PER_GEV, self.pointsy / ERGS_PER_GEV, color=color,
            label=self.get_name(), marker=shape, facecolors=fcolor)
        thisx = []
        thisy = []
        up_err = []
        down_err = []
        lolims=[]
        uplims=[]
        for x in self.pointsx:
            y = self.get_point_log(x)
            u = self.get_up_bar_log(x)
            d = self.get_down_bar_log(x)
            thisx.append(x)
            thisy.append(y)
            if u is None and d is None:
                up_err.append(0)
                down_err.append(0)
                uplims.append(False)
                lolims.append(False)
            elif u is None:
                if d < Y_LIM:
                    down_err.append(y*ERROR_BAR_SIZE)
                    lolims.append(True)
                else:
                    down_err.append(y-d)
                    lolims.append(False)
                up_err.append(y*ERROR_BAR_SIZE)
                uplims.append(True)
            elif d is None:
                up_err.append(u-y)
                uplims.append(False)
                down_err.append(y*ERROR_BAR_SIZE)
                lolims.append(True)
            else:
                up_err.append(u-y)
                uplims.append(False)
                if d < Y_LIM:
                    down_err.append(y*ERROR_BAR_SIZE)
                    lolims.append(True)
                else:
                    down_err.append(y-d)
                    lolims.append(False)
        thisx = np.asarray(thisx)
        thisy = np.asarray(thisy)
        down_err = np.asarray(down_err)
        up_err = np.asarray(up_err)
        ax.errorbar(10**thisx / ERGS_PER_GEV, thisy / ERGS_PER_GEV, yerr=[down_err / ERGS_PER_GEV, up_err / ERGS_PER_GEV], color=color,
            linestyle='none', elinewidth=1, lolims=uplims, uplims=lolims)
        ax.set_title(self.get_name())
        ax.set_xlim(0.1, 100)
        ax.set_yscale("log", nonpositive='clip')
        ax.set_xscale("log", nonpositive='clip')

    def display_calore(self, ax, l_min, l_max, show_all=False):
        self._fit_calore(l_min, l_max, False)

        if not show_all:
            calorex = np.linspace(np.log10(l_min * ERGS_PER_GEV),
                np.log10(l_max * ERGS_PER_GEV), 100)
        else:
            calorex = np.linspace(np.min(self.pointsx),
                np.max(self.pointsx), 100)
        calorey = calore_power_law(calorex, self.calore_norm)

        ax.plot(10**calorex / ERGS_PER_GEV, calorey / ERGS_PER_GEV, color='C1', linestyle='--')

    def display_power_law(self, ax, l_min, l_max, show_all=False):
        self._fit_power_law(l_min, l_max, False)

        if not show_all:
            fitx = np.linspace(np.log10(l_min * ERGS_PER_GEV),
                np.log10(l_max * ERGS_PER_GEV), 100)
        else:
            fitx = np.linspace(np.min(self.pointsx),
                np.max(self.pointsx), 100)
        fity = power_law(fitx, self.fit_norm, self.fit_alpha_below,
            self.fit_alpha_above, self.fit_l_break) / LM_FIT_SCALE

        ax.plot(10**fitx / ERGS_PER_GEV, fity / ERGS_PER_GEV, color='C2', linestyle='-.')

    def get_x_label(self):
        return "$E_\\gamma$ [GeV]"

    def get_y_label(self):
        return "$F_\\gamma$ [GeV / cm$^2$ / s]"
