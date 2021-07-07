DI_MAURO = 0
FERMILAB = 1
CALORE = 2
TEST = 3
ERG_PER_GEV = 0.00160218
FLUX_TO_LUM = 1.1093417307914119e-46
NPTF_FLUX = [1.76e-10, 6.80e-9]

TYPE = DI_MAURO

if TYPE == FERMILAB:
    L_BREAK = 0.00516323313751923
    ALPHA_ABOVE = 2.7922607115665783
    ALPHA_BELOW = 1.561888668638258
if TYPE == CALORE:
    L_BREAK = 2.06 * ERG_PER_GEV
    ALPHA_ABOVE = 2.63
    ALPHA_BELOW = 1.42
if TYPE == DI_MAURO:
    L_BREAK = 0.0020217250786090323
    ALPHA_ABOVE = 2.565547079994774
    ALPHA_BELOW = 1.0337152835197816
if TYPE == TEST:
    L_BREAK = 0.1
    ALPHA_ABOVE = 1.2
    ALPHA_BELOW = 1.4

def integrate(start, stop):
    topAlpha = ALPHA_ABOVE if stop > L_BREAK else ALPHA_BELOW
    bottomAlpha = ALPHA_BELOW if start < L_BREAK else ALPHA_ABOVE
    return L_BREAK / (1 - bottomAlpha) * (1 - (start / L_BREAK)**(1 - bottomAlpha)) + \
        L_BREAK / (1 - topAlpha) * ((stop / L_BREAK)**(1 - topAlpha) - 1) # Math confirmed 7 Jul 2021


def lintegrate(start, stop):
    topAlpha = ALPHA_ABOVE if stop > L_BREAK else ALPHA_BELOW
    bottomAlpha = ALPHA_BELOW if start < L_BREAK else ALPHA_ABOVE
    return L_BREAK**2 / (2 - bottomAlpha) * (1 - (start / L_BREAK)**(2 - bottomAlpha)) + \
        L_BREAK**2 / (2 - topAlpha) * ((stop / L_BREAK)**(2 - topAlpha) - 1) # Math confirmed 7 Jul 2021

def getPhotonEnergy(start, stop):
    return lintegrate(start, stop) / integrate(start, stop) # Units of ergs

if TYPE == FERMILAB:
    spectrumName = "Fermilab"
if TYPE == CALORE:
    spectrumName = "Calore"
if TYPE == DI_MAURO:
    spectrumName = "Di Mauro"
if TYPE == TEST:
    spectrumName = "test"
print("Using {} spectrum".format(spectrumName))

print("Energy per photon for 0.1-100 GeV spectrum (erg):", getPhotonEnergy(0.1 * ERG_PER_GEV, 100 * ERG_PER_GEV))
print("Energy per photon for NPTF spectrum (erg):", getPhotonEnergy(1.893 * ERG_PER_GEV, 11.943 * ERG_PER_GEV))
print()
print("NPTF lum for 0.1-100 GeV (erg/s):", [f * getPhotonEnergy(0.1 * ERG_PER_GEV, 100 * ERG_PER_GEV)/FLUX_TO_LUM for f in NPTF_FLUX])
print("NPTF lum for NPTF spectrum (erg/s):", [f * getPhotonEnergy(1.893 * ERG_PER_GEV, 11.943 * ERG_PER_GEV)/FLUX_TO_LUM for f in NPTF_FLUX])
