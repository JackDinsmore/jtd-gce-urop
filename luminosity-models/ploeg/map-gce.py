from math import exp, log10
from matplotlib import pyplot as plt

POINTS = 200
X_RANGE = [0.275, 51.9] # GeV
X = []
Y = []

K = 1 # Normalization constant
PI = 3.1415926535
GAMMA = 0.163
E_CUT = 10**0.039 

def dNdE(E):
    return K * E**(-GAMMA) * exp(-(E / E_CUT)**(2/3))

def EsquareddNdE(E):
    return E * E * dNdE(E)

def populate():
    x = X_RANGE[0]
    while x < X_RANGE[1]:
        X.append(x)
        Y.append(EsquareddNdE(x))
        x += (X_RANGE[1] - X_RANGE[0]) / POINTS

def getGCEFlux():
    E = 0.1
    integral = 0
    dE = 0.01
    while E < 100:
        integral += E * dNdE(E) * dE
        E += dE
    return integral

print("Predicted flux of GCE: {0} ergs/sec".format(getGCEFlux()))

populate()
plt.xscale('log')
plt.yscale('log')
plt.title("GCE flux distribution predicted by Ploeg")
plt.xlabel("Energy")
plt.ylabel("E^2 dN / dE (arbitrary units for now)")
plt.plot(X, Y)
plt.show()