from math import sqrt, exp, pi

IN_FILE = "disk.csv"
OUT_FILE = "power-law.csv"
LOG_L0 = 32.206
SIGMA = 0.70585

MIN_LOG_L = 28
MAX_LOG_L = 37
SKIP = 0.001

ALPHA = 1.93
L_MAX = 1e35

def getPowerLaw(logl):
    return 1e55 * (10**logl)**(-ALPHA) * exp(-(10**logl) / L_MAX)

logl = MIN_LOG_L
newData = ''
while logl < MAX_LOG_L:
    newData += str(logl) + "," + str(getPowerLaw(logl)) + "\n"
    logl += SKIP

f = open(OUT_FILE, 'w')
f.write(newData)
f.close()