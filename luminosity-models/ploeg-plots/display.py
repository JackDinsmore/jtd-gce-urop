from matplotlib import pyplot as plt
from math import log10

DISK=0
BOXY_BULGE = 1
NUCLEAR_BULGE = 2

INTEGRATION_LOG_STEP = 0.001

class LuminosityFunction:
    def __init__(self, popType):
        self.popType = popType
        if popType not in [0, 1, 2]:
            raise Exception("Population type {0} is not available. Options are DISK, BOXY_BULGE, and NUCLEAR_BULGE.".format(popType))
        self.data =[]
        f = open(self.getFile())
        lines = f.read().split('\n')
        f.close()
        for line in lines:
            if line == '': continue
            x, y = line.split(",")
            self.data.append([float(x), float(y)])
        
        self.leftPoint = self.data[0]
        self.rightPoint = self.data[-1]

        self.base = None
        if popType != DISK:
            self.base = LuminosityFunction(DISK)

    def __call__(self, l):
        if self.base == None:
            return self.getThisLogValue(log10(l))
        else:
            return self.getThisLogValue(log10(l)) + self.base.getThisLogValue(log10(l))

    def getThisLogValue(self, l):
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
        names = ["Disk", "Boxy bulge", "Nuclear bulge"]
        return names[self.popType]

    def getFile(self):
        prepend = "C:/Users/goods/Dropbox (MIT)/GCE UROP/luminosity-models/ploeg-plots/"
        files = ["disk.csv", "boxy-bulge.csv", "nuclear-bulge.csv"]
        return prepend + files[self.popType]

    def display(self):
        x = []
        y = []
        step = 0.1
        logl = self.leftPoint[0]
        while logl < self.rightPoint[0]:
            x.append(10**logl)
            y.append(self.__call__(10**logl))
            logl += step
        plt.title("{0} luminosity function".format(self.getName()))
        plt.xlabel("Luminosity [erg s^(-1)]")
        plt.ylabel("Probability density")
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
            integral += self.__call__(10**logl) * (10**(logl + INTEGRATION_LOG_STEP) - 10**(logl))
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
            integral += 10**logl * self.__call__(10**logl) * (10**(logl + INTEGRATION_LOG_STEP) - 10**(logl))
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
    plotAll()