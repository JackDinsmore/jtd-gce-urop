# Refining the Millisecond Pulsar Model for the Galactic Center Excess

## Description

This is a repository of plots and code made by Jack Dinsmore, started in October of 2020. Please see _proposal/proposal.pdf_ for an introduction to the project.

The _summaries/_ directory contains roughly monthly-made reports on the progress I've made and the purpose of the plots I produced.

The _luminosity\_models/_ directory contains the code itself, designed to analyze the property of the GCE as modeled by the following three luminosity functions:

* **Power law distribution.** The simplest luminosity function. This was used by the Fermilab paper Zhong et al., which is cited in the project proposal.

* **Log normal distribution.** This luminosity function is slightly more complex and seems to fit the data better. It was used by Hooper et al., cited in the October 2020 summary

* **Ploeg distribution.** In Ploeg et al., cited in the October 2020 summary, a detailed analysis of neutron stars is used to generate a numerical luminosity function which has been manually extracted from that paper. It very closely resembles a log-normal distribution.

In addition to these three luminosity functions, used to model the GCE population, we also need to model the _Fermi_ telescope, whose observations are a central feature of this project. We use the following two models:

* **Step function sensitivity.** We assume that _Fermi_ sees all the stars below a certain threshold (10<sup>34</sup> ergs per second)

* **Error function sensitivity.** We use a more detailed model of the sensitivity which closely resembles a step function, taken from Ploeg et al. equation 2.34. 
<div align="center"> <img src="https://render.githubusercontent.com/render/math?math=P(F_{th} \leq F)=\frac{1}{2} \left[1 %2B \erf \left( \frac{\log_{10}(F) - (\log_{10}(\mu_{th}(l, b)) %2B K_{th}}{\sqrt{2}\sigma_{th}}\right)\right]"></div>
Foo