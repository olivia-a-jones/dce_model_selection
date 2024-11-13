# dce_model_selection

![dce models](images/Figure1.png)

Scripts used to perform tracer kinetic model selection on simulated and in-vivo DCE-MRI data.
All quantities, processes and model definitions are [OSIPI CAPLEX compliant](https://doi.org/10.1002/mrm.29840).<sup>1</sup> CAPLEX definitions can be accessed by clicking on quantity, process or model hyperlinks.

---
### 1. Required software
Madym<sup>2</sup> and Madym's Python wrappers are required to run these scripts and can be downloaded [here](https://gitlab.com/manchester_qbi/manchester_qbi_public/madym_cxx).

--- 
### 2. Simulating DCE-MRI signal-time curves using an existing VIF and performing model selection. 
This script does the following:
- Simulates signal time-series with user-specified tissue parameters with the [2CXM](https://osipi.github.io/OSIPI_CAPLEX/perfusionModels/#2CXM) using our smoothed group-averaged VIF ([the indicator concentration time series for blood plasma](https://osipi.github.io/OSIPI_CAPLEX/quantities/#C)) from control participants. The ground-truth PS and the amount of noise added to the curve are varied to produce a grid of 20 x 20 noisy time series. 
  - The group-averaged VIF is provided in group_averaged_VIF.txt
  - Details of the participants and DCE-MRI acquisition relating to the VIF can be found in [this paper](https://doi.org/10.3389/fphys.2020.593026)).
- Fits an [Extended Tofts model](https://osipi.github.io/OSIPI_CAPLEX/perfusionModels/#ETM) of indicator exchange, a [Patlak model](https://osipi.github.io/OSIPI_CAPLEX/perfusionModels/#Patlak) of indicator uptake, and an intravascular model.
- Selects the best fitting model for each time-series in the grid using the [Akaike Information Criterion](https://osipi.github.io/OSIPI_CAPLEX/quantities/#AIC).
- Saves the best-fitting model grid, the fitted K<sup>trans</sup> grid, 

---
### 3. Performing model selection on in-vivo DCE-MRI data. 
This script is designed to process our scans (details of which can be found in the supplementary materials of [this paper](https://doi.org/10.3389/fphys.2020.593026)), and may need to be edited to be run on other DCE-MRI datasets. An example dataset can be requested by [email](olivia.jones-4@manchester.ac.uk).

---
### References
1. Dickie BR, Ahmed Z, Arvidsson J, et al. A community-endorsed open-source lexicon for contrast agent-based perfusion MRI: A consensus guidelines report from the ISMRM Open Science Initiative for Perfusion Imaging (OSIPI). Magn Reson Med. Published online October 13, 2023. doi:10.1002/mrm.29840
2. Berks M, Parker GJM, Little R, Cheung S. Madym: A C++ toolkit for quantitative DCE-MRI analysis. Published online 2021. https://doi.org/10.5281/zenodo.5176079
