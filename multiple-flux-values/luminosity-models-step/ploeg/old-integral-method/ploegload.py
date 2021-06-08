from matplotlib import pyplot as plt
from math import log10, exp

DEBUG = True

DISK=0
BOXY_BULGE = 1
NUCLEAR_BULGE = 2
LOG_NORMAL_FIT = 3

INTEGRATION_LOG_STEP = 0.001

PREPEND = "C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models-step/ploeg/data/" # Change if you use this on a different system

class LuminosityFunction:
    def __init__(self, popType):
        self.popType = popType
        if popType not in [0, 1, 2] and not DEBUG:
            raise Exception("Population type {0} is not available. Options are DISK, BOXY_BULGE, and NUCLEAR_BULGE.".format(popType))
        self.data =[]
        f = open(self.getFile())
        lines = f.read().split('\n')
        f.close()
        for line in lines:
            if line == '': continue
            x, y = line.split(",")# log L, dn/dlog L
            newY = log10(exp(1)) / 10**float(x) * float(y) #10**float(x)*float(y)*log(10)
            self.data.append([float(x), newY])# log L, dn/ dL
        
        self.leftPoint = self.data[0]
        self.rightPoint = self.data[-1]

        self.base = None
        if popType not in [DISK, LOG_NORMAL_FIT]:
            self.base = LuminosityFunction(DISK)

    def __call__(self, l): # L given in decimal
        return self.callWithLogLum(log10(l))

    def callWithLogLum(self, logl):
        if self.base == None:
            return self.getThisLogValue(logl)
        else:
            return self.getThisLogValue(logl) + self.base.getThisLogValue(logl)


    def getThisLogValue(self, l): # Log10 l given.
        if l <= self.leftPoint[0]:
            return self.leftPoint[1]
        if l >= self.rightPoint[0]:
            return self.rightPoint[1]
        
        # Binary search for l:
        leftIndex = 0
        rightIndex = len(self.data) - 1
        while True:
            midIndex = (rightIndex + leftIndex)//2
            if self.data[midIndex][0] > l:
                rightIndex = midIndex
            elif self.data[midIndex][0] < l:
                leftIndex = midIndex
            else:
                return self.data[midIndex][1]
            if rightIndex - leftIndex <= 1:
                assert leftIndex != rightIndex
                return self.data[leftIndex][1] + (self.data[rightIndex][1] - self.data[leftIndex][1]) \
                    * (l - self.data[leftIndex][0]) \
                        / (self.data[rightIndex][0] - self.data[leftIndex][0])

    def getName(self):
        names = ["Disk", "Boxy bulge", "Nuclear bulge", "Log normal fit"]
        return names[self.popType]

    def getFile(self):
        files = ["disk.csv", "boxy-bulge.csv", "nuclear-bulge.csv", "log-normal-fit.csv"]
        return PREPEND + files[self.popType]

    def display(self):
        x = []
        y = []
        step = 0.1
        logl = self.leftPoint[0]
        while logl < self.rightPoint[0]:
            x.append(10**logl)
            y.append(self.callWithLogLum(logl))
            logl += step
        plt.title("{0} luminosity function".format(self.getName()))
        plt.xlabel("Luminosity [erg s$^{-1}$]")
        plt.ylabel("dN/dL")
        plt.xscale('log')
        plt.plot(x, y, label=self.getName())

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

    def lintegrate(self, minL=None, maxL=None):
        if minL == None:
            minL = 10**self.data[0][0]
        if maxL == None:
            maxL = 10**self.data[-1][0]
        if maxL == minL:
            return 0
        if maxL < minL:
            return -self.lintegrate(minL=maxL, maxL=minL)

        logl = log10(minL)
        integral = 0
        while logl < log10(maxL):
            # Trapezoidal integration
            xLow = 10**(logl)
            xHigh = 10**(logl + INTEGRATION_LOG_STEP)
            yLow = self.callWithLogLum(logl)
            yHigh = self.callWithLogLum(logl + INTEGRATION_LOG_STEP)
            integral += (xLow * yLow + xHigh * yHigh) / 2 * (xHigh - xLow)
            logl += INTEGRATION_LOG_STEP
        return integral


def plotAll():
    plt.title("All luminosity functions")
    plt.xlabel("Luminosity [erg s^(-1)]")
    plt.ylabel("Probability density")
    plt.xscale('log')

    functions = [LuminosityFunction(DISK),
        LuminosityFunction(BOXY_BULGE),
        LuminosityFunction(NUCLEAR_BULGE)]
    for f in functions:
        x = []
        y = []
        step = 0.1
        logl = f.leftPoint[0]
        while logl < f.rightPoint[0]:
            x.append(10**logl)
            y.append(f(10**logl))
            logl += step
        plt.plot(x, y, label=f.getName())
    plt.legend()
    plt.show()
    plt.savefig("all-plots.png")

if __name__ == "__main__":
    f = LuminosityFunction(LOG_NORMAL_FIT)
    norm = f.integrate()
    print(f.integrate(10**31.28))
    print(f.lintegrate(10**31.28))
    print(f.integrate(10**31.28) / norm)
    print(f.lintegrate(10**31.28) / norm)