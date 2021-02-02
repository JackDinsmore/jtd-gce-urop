# Refining the Millisecond Pulsar Model for the Galactic Center Excess

## Description

This is a repository of plots and code made by Jack Dinsmore, started in October of 2020. Please see _proposal/proposal.pdf_ for an introduction to the project.

## Directories

The _compare-spectra/_ directory calculates the total luminosity of the GCE as a function of the radiation energies included. Common ranges include 0.1-100 GeV, yielding a total luminosity of 6.37e34 ergs per second, but other ranges can be accounted for with this code.

The _flux-histograms/_ directory displays a histogram of the luminosities and fluxes observed from MSPs in the GCE.

The _luminosity-models*/_ directory contains the code itself, designed to analyze the property of the GCE as modeled by the following luminosity functions:

* **Power law distribution.** The simplest luminosity function. This was used by the Fermilab paper Zhong et al., which is cited in the project proposal.

* **Log normal distribution.** This luminosity function is slightly more complex and seems to fit the data better. It was used by Hooper et al., cited in the October 2020 summary

* **Ploeg distribution.** In Ploeg et al., cited in the October 2020 summary, a detailed analysis of neutron stars is used to generate a numerical luminosity function which has been manually extracted from that paper. It very closely resembles a log-normal distribution.

* **NPTF distribution.** A broken power law luminosity function whose parameters are fixed by the Non-Poissonian Template Fitting method (NPTF) outlined in S K Lee et al. (2016).

In addition to these three luminosity functions, used to model the GCE population, we also need to model the _Fermi_ LAT, whose observations are a central feature of this project. We use the following two models:

* **Step function sensitivity.** We assume that the LAT sees all the stars with luminosity above 10<sup>34</sup> ergs per second. This is carried out in the _luminosity-models-step/_ directory.

* **Numerical sensitivity.** Using a map of the sensitivity of the LAT at all points in the sky, we assuming that the point sources assumed to make up the GCE are distributed according to a Navarro-Frank-White (NFW) squared profile and calculate the number and fraction of total luminosity observed in each pixel of the sensitivity map, adding together the results to achieve a more accurate result than the step function sensitivity model. This is carried out in the _luminosity-models-position/_ directory, and the sensitivity map is contained in _sensitivity_.

The _summaries/_ directory contains roughly monthly-made reports on the progress I've made and the purpose of the plots I produced.
