import numpy as np

ALPHA_BELOW = 1.1
ALPHA_ABOVE = 2.5
SCALE = 1e-7;
ERGS_PER_GEV = 0.00160218

L_BREAK=1.23

E_MIN = np.log10(0.1)
E_MAX = np.log10(100)
N_POINTS = 1000


def power_law(x, scale, alpha_below, alpha_above, l_break):
    alpha = np.full_like(x, alpha_above)
    alpha[10**x < l_break] = alpha_below
    return (10**x)**2 * scale * (10**x / l_break)**(-alpha)

f = open("points.csv", 'w')
for logx in np.linspace(E_MIN, E_MAX, N_POINTS):
    f.write(str(logx) + ", " + str(np.log10(power_law(logx, SCALE, ALPHA_BELOW, ALPHA_ABOVE, L_BREAK)))+"\n")
f.close()

f = open("down-bars.csv", 'w')
for logx in np.linspace(E_MIN, E_MAX, N_POINTS):
    f.write(str(logx) + ", " + str(np.log10(power_law(logx, SCALE, ALPHA_BELOW, ALPHA_ABOVE, L_BREAK))*0.99)+"\n")
f.close()

f = open("up-bars.csv", 'w')
for logx in np.linspace(E_MIN, E_MAX, N_POINTS):
    f.write(str(logx) + ", " + str(np.log10(power_law(logx, SCALE, ALPHA_BELOW, ALPHA_ABOVE, L_BREAK))*1.01)+"\n")
f.close()
