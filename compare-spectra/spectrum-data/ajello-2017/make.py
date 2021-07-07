import numpy as np
import matplotlib.pyplot as plt

EPSILON = 1e-15

def get_arrays(filename):
    f = open(filename, 'r')
    xs = []
    ys = []
    for line in f.readlines():
        if line == '': continue
        x, y = line.split(", ")
        xs.append(float(x))
        ys.append(float(y))
    return np.asarray(xs), np.asarray(ys)

px, py = get_arrays("raw-points.csv")
_, bary = get_arrays("raw-bars.csv")

plt.errorbar(px, py, yerr=bary)
plt.xscale("log")
plt.yscale("log")
plt.ylim(3e-9, 2e-6)

points = open("points.csv", 'w')
up_bars = open("up-bars.csv", 'w')
down_bars = open("down-bars.csv", 'w')
for i in range(len(px)):
    if (py[i] < 0): continue
    points.write(str(np.log10(px[i]))+", "+str(np.log10(py[i]))+'\n')
    up_bars.write(str(np.log10(px[i]))+", "+str(np.log10(py[i]+bary[i]))+'\n')
    down_bars.write(str(np.log10(px[i]))+", "+str(np.log10(max(py[i]-bary[i], EPSILON)))+'\n')

points.close()
up_bars.close()
down_bars.close()

plt.show()
