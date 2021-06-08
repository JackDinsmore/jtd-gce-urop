OLD = 0
FERMILAB = 1
CALORE = 2
ERG_PER_GEV = 0.00160218
SIZE_OF_ROI = 0.42882 # steradians
L_BREAK = 2.06 * ERG_PER_GEV
ALPHA_ABOVE = 2.63
ALPHA_BELOW = 1.42

TYPE = FERMILAB

# The type doesn't actually matter; it doesn't come into the calcualtion

SCALE = 0 # Units in erg / cm^2 / s
if TYPE == OLD:
    raise Exception("Old values unimplemented")
if TYPE == FERMILAB:
    SCALE = 25.87693692 * ERG_PER_GEV * SIZE_OF_ROI
if TYPE == CALORE:
    SCALE = 44.78745446 * ERG_PER_GEV * SIZE_OF_ROI

def integrate(start, stop):
    topAlpha = ALPHA_ABOVE if stop > L_BREAK else ALPHA_BELOW
    bottomAlpha = ALPHA_BELOW if start < L_BREAK else ALPHA_ABOVE
    return L_BREAK / (1 - bottomAlpha) * (1 - (start / L_BREAK)**(1 - bottomAlpha)) + \
        L_BREAK / (1 - topAlpha) * ((stop / L_BREAK)**(1 - topAlpha) - 1)
    

def lintegrate(start, stop):
    topAlpha = ALPHA_ABOVE if stop > L_BREAK else ALPHA_BELOW
    bottomAlpha = ALPHA_BELOW if start < L_BREAK else ALPHA_ABOVE
    return L_BREAK**2 / (2 - bottomAlpha) * (1 - (start / L_BREAK)**(2 - bottomAlpha)) + \
        L_BREAK**2 / (2 - topAlpha) * ((stop / L_BREAK)**(2 - topAlpha) - 1)

def getPhotonEnergy(start, stop):
    return lintegrate(start, stop) / integrate(start, stop) # Units of ergs

spectrumName = "old"
if TYPE == FERMILAB:
    spectrumName = "Fermilab"
if TYPE == CALORE:
    spectrumName = "Calore"
print("Using {} spectrum".format(spectrumName))

print("Energy per photon for 0.1-100 GeV spectrum:", getPhotonEnergy(0.1 * ERG_PER_GEV, 100 * ERG_PER_GEV))
print("Energy per photon for NPTF spectrum:", getPhotonEnergy(1.893 * ERG_PER_GEV, 11.943 * ERG_PER_GEV))