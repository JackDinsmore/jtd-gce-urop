from matplotlib import pyplot as plt
from math import log10, exp

DISK=0
LOG_NORMAL_FIT = 3

INTEGRATION_LOG_STEP = 0.001

PREPEND = "C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models-step/ploeg/data/" # Change if you use this on a different system

class LuminosityFunction:
    def __init__(self, popType):
        self.popType = popType
        if popType not in [0, 3]:
            raise Exception("Population type {0} is not available. Options are DISK and LOG_NORMAL_FIT.".format(popType))
        self.xs =[]
        self.ys =[]
        self.logxs =[]
        self.logys =[]
        f = open(self.getFile())
        lines = f.read().split('\n')
        f.close()
        print(len(lines))
        for line in lines:
            if line == '': continue
            x, y = line.split(",")# log L, dn/dlog L
            newY = log10(exp(1)) / 10**float(x) * float(y) #10**float(x)*float(y)*log(10)
            self.xs.append(10**float(x))
            self.ys.append(newY)
            self.logxs.append(float(x))
            self.logys.append(float(y))

    def getName(self):
        names = ["Disk", "Boxy bulge", "Nuclear bulge", "Log normal fit"]
        return names[self.popType]

    def getFile(self):
        files = ["disk.csv", "boxy-bulge.csv", "nuclear-bulge.csv", "log-normal-fit.csv"]
        return PREPEND + files[self.popType]

    def display(self):
        step = 0.1
        plt.title("{0} luminosity function".format(self.getName()))
        plt.xlabel("Luminosity [erg s$^{-1}$]")
        plt.ylabel("dN/dL")
        plt.xscale("log")
        plt.plot(self.xs, self.logys, label=self.getName())

    def integrate(self, minL=None, maxL=None):
        if minL == None:
            minL = 10**self.data[0][0]
        if maxL == None:
            maxL = 10**self.data[-1][0]
        if maxL == minL:
            return 0
        if maxL < minL:
            return -self.integrate(minL=maxL, maxL=minL)

        logl = log10(minL)
        integral = 0
        while logl < log10(maxL):
            # Trapezoidal integration
            xLow = 10**(logl)
            xHigh = 10**(logl + INTEGRATION_LOG_STEP)
            yLow = self.callWithLogLum(logl)
            yHigh = self.callWithLogLum(logl + INTEGRATION_LOG_STEP)
            integral += (yLow + yHigh) / 2 * (xHigh - xLow)
            logl += INTEGRATION_LOG_STEP
        return integral



    def integrate(self, start):
        i = 0
        while i < len(self.xs) and self.xs[i] < start:
            i+= 1
        assert(i < len(self.xs))
        assert(i > 0)
        frac = (start - self.xs[i - 1]) / (self.xs[i] - self.xs[i - 1])
        startY = self.ys[i - 1] + (self.ys[i] - self.ys[i - 1]) * frac
        integral = 0.5 * (startY + self.ys[i]) * (self.xs[i] - start)
        while i < len(self.xs) - 1:
            integral += 0.5 * (self.ys[i] + self.ys[i + 1]) * (self.xs[i + 1] - self.xs[i])
            i+= 1
        return integral

    def lintegrate(self, start):
        i = 0
        while i < len(self.xs) and self.xs[i] < start:
            i+= 1
        assert(i < len(self.xs))
        assert(i > 0)
        frac = (start - self.xs[i - 1]) / (self.xs[i] - self.xs[i - 1])
        startY = self.ys[i - 1] + (self.ys[i] - self.ys[i - 1]) * frac
        integral = 1 / 6.0 * (-start * start * (self.ys[i] + 2 * startY) + start * self.xs[i] * (startY - self.ys[i]) + self.xs[i] * self.xs[i] * (2 * self.ys[i] + startY))
        while i < len(self.xs) - 1:
            integral += 1 / 6.0 * (self.ys[i] * (-2 * self.xs[i] * self.xs[i] + self.xs[i] * self.xs[i + 1] + self.xs[i + 1] * self.xs[i + 1])
				+ self.ys[i + 1] * (-self.xs[i] * self.xs[i] - self.xs[i] * self.xs[i + 1] + 2 * self.xs[i + 1] * self.xs[i + 1]))
            i+= 1   
        return integral



if __name__ == "__main__":
    f = LuminosityFunction(LOG_NORMAL_FIT)
    norm = f.integrate()
    print(f.integrate(10**31.28))
    print(f.lintegrate(10**31.28))
    print(f.integrate(10**31.28) / norm)
    print(f.lintegrate(10**31.28) / norm)