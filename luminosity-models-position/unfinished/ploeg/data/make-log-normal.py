from math import sqrt, exp, pi

IN_FILE = "disk.csv"
OUT_FILE = "log-normal-fit.csv"
LOG_L0 = 32.206
SIGMA = 0.70585

MIN_LOG_L = 27
MAX_LOG_L = 36
SKIP = 0.001

ALPHA = 1.93
L_MAX = 1e35

def getLogNormal(logl):
    return 1 / (SIGMA * sqrt(2 * pi)) * exp(-(logl - LOG_L0)**2 / (2 * SIGMA**2))

logl =MIN_LOG_L
newData = ''
while logl < MAX_LOG_L:
    newData += str(logl) + "," + str(getLogNormal(logl)) + "\n"
    logl += SKIP

f = open(OUT_FILE, 'w')
f.write(newData)
f.close()