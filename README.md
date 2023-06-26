Structural and dynamical quality assessment of gap-filled sea surface temperature products
============

Authors
--------
[Cristina González-Haro](https://crocha700.github.io)<sup>1</sup>, [Jordi Isern-Fontanet](http://www-pord.ucsd.edu/~sgille/)<sup>1</sup>
, [Antonio Turiel](http://tryfan.ucsd.edu)<sup>1</sup>,
 [Christopher J. Merchant](https://science.jpl.nasa.gov/people/Menemenlis/)<sup>2</sup> and  [Peter Cornillon](https://science.jpl.nasa.gov/people/Menemenlis/)<sup>3</sup> .



1: [Institut de Cièncias del Mar (ICM-CSIC)](https://icm.csic.es),
[Barcelona Expert Center on Remote Sensing](https://bec.icm.csic.es), Barcelona, Spain.

2: [National Centre for Earth Observation and Department of Meteorology](https://www.nceo.ac.uk/), [University of Reading](https://www.reading.ac.uk/), Reading, England, UK

3: Graduate School of Oceanography, [University of Rhode Island](https://www.uri.edu/), Narragansett, RI, USA.


Key Points
----------
  -  We proposed a new diagnostic to assess the structural and dynamical properties of sea surface temperature ( SST ) products
  - This diagnostic is based on multifractal theory of turbulence and consists in computing the singularity spectrum.
  - The different schemes used to produce gap-filled SST products may contribute to the loss of dynamical information or structural coherence.

Abstract
--------
Previous studies that intercompared global L4 SST analyses were centered on the assessment of the accuracy and bias of SST by comparing them with independent near-surface Argo profile temperature data. This type of assessment is centered on the absolute value of SST rather than on SST spatial properties (gradients), which is more relevant to the study of oceanographic features (e.g., fronts, eddies, etc) and ocean dynamics. Here,  we  use, for the first time, the spectrum of singularity exponents to assess the structural and dynamical quality of different L4 gap-filled products based on the multifractal theory of turbulence. Singularity exponents represent the geometrical projection of the turbulence cascade, and its singular spectrum can be related to the PDF of the singularity exponents normalized by the scale. Our results reveal that the different schemes used to produce the L4 SST products generate different singularity spectra, which are then used to identify a potential loss of dynamical information or structural coherence. This new diagnostic constitutes a valuable tool to assess the structural quality of SST products and can support data satellite SST producers efforts to improve the interpolation schemes used to generate gap-filled SST products. 

Status
----------
  The paper is submitted to [Earth Space Sicence](https://agupubs.onlinelibrary.wiley.com/journal/23335084). Comments, questions, and suggestions are extremely welcome  and warmly appreciated. Feedback can be submitted via e-mail to  Cristina González-Haro (cgharo@icm.csic.es).

Code
----
The analysis for this paper has been performed on ICM's computing cluster, *Gaia*. We first compute the singularity exponents of 5 different SST products using the algorithm proposed by [Pont et al. (2013)](https://www.tandfonline.com/doi/abs/10.1080/00207160.2012.748895). Once the computations are performed on *Gaia*, the output files are saved in netCDF4 and  are publicably available at [DigitalCSIC](https://digital.csic.es/handle/10261/309548).
two small pieces of code developed by the first author and available on github:  [llctools](https://github.com/crocha700/llctools) and [pyspec](https://github.com/pyspec/pyspec). Those codes leverage on the [Scientific
ython stack](https://www.scipy.org/install.html). Specific processing and plotting codes are available on Jupyter [notebooks](https://github.com/crocha700/UpperOceanSeasonality/blob/master/notebooks/index.ipynb).



Data
------
The input data to produce the singularity analysis described in the paper are the singularity exponents publicaly available at [DigitalCSIC](https://digital.csic.es/handle/10261/309548). The source SST products are available at:
  - [SST_AMSRE2_REMSS](https://data.remss.com/amsr2/ocean/L3/)
  - [SST_CMC](https://podaac.jpl.nasa.gov/dataset/CMC0.1deg-CMC-L4-GLOB-v3.0)
  - [OSTIA](https://resources.marine.copernicus.eu/product-detail/SST_GLO_SST_L4_NRT_OBSERVATIONS_010_001/INFORMATION)
  - [SST_CCI](https://resources.marine.copernicus.eu/product-detail/SST_GLO_SST_L4_REP_OBSERVATIONS_010_024)
  - [SST_MUR](https://podaac.jpl.nasa.gov/dataset/MUR-JPL-L4-GLOB-v4.1)

Support
-------
This work was supported by the DEMON project (PID2021-123457OB-C21) and the INTERACT project (PID2020-114623RB-C31) funded by the Spanish Ministerio de Ciencia e Innovación; the Agencia Estatal de Investigación (10.13039/501100011033); and the European Regional Development Fund, EU. This work also acknowledges the ‘‘Severo Ochoa Centre of Excellence’’ accreditation (CEX2019-000928-S). This work is moreover a contribution to [CSIC PTI Teledetect](https://pti-teledetect.csic.es/).

DOI
-------
[![DOI](https://zenodo.org/badge/657566970.svg)](https://zenodo.org/badge/latestdoi/657566970)
