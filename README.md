[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17962642.svg)](https://doi.org/10.5281/zenodo.17962642)
# The topology of software risk in scientific models

[Arnald Puy](https://www.arnaldpuy.com/), Federico Díaz, Ulrike Proske, Olivia Richards,
Seth N. Linga, Samuel Flinders, Carmen Aguiló-Rivera, and Fernando G. Tinetti.

This study proposes a framework to identify risky paths in scientific models; that is, 
sequences of function calls whose potential failure can have a larger cascading effect into
other parts of the software. We illustrate our method with a synthetic example and by applying it
to fourteen global land use and hydrological models, listed in the Models section below.

## Abstract

*Scientific models are becoming complex software systems whose structural fragility 
remains largely hidden. To address this issue we present a domain-agnostic framework 
that integrates software metrics with network analysis, representing each model’s 
source code as a directed call graph where nodes and edges correspond to functions 
and their calls. By combining cyclomatic complexity, call centrality and connectivity 
metrics we compute node- and path-level risk scores that quantify how software faults 
may propagate through the code architecture. This approach exposes "risk highways", 
chains of interdependent functions where local faults are most likely to cascade. 
We illustrate our approach on fourteen large-scale hydrological and land-surface models. 
Our framework reveals the hidden topology of software risk and provides a quantitative 
basis for assessing, comparing and improving the reliability of scientific models.*

## Models

The models used in this study are the following:

* [CTSM](https://github.com/ESCOMP/CTSM)  - Community Terrestrial Systems Model.    
* [CWatM](https://github.com/iiasa/CWatM) - Community Water Model.      
* [DBH](https://hydro.iis.u-tokyo.ac.jp/DBH/index_files/Page394.htm) - Distributed Biosphere-Hydrological Model.
* [GR4J](https://github.com/EdgarEspitia/GR4J) - Génie Rural à 4 paramètres Journalier.       
* [H08](https://github.com/h08model/H08) - H08 Global Hydrological Model.
* [HBV](https://github.com/johnrobertcraven/hbv_hydromodel) - Hydrologiska Byråns Vattenbalansavdelning model.      
* HydroPy - (Hydrological model implemented in Python) 
* [HYPE](https://sourceforge.net/projects/hype/files/) - Hydrological Predictions for the Environment.
* [MHM](https://zenodo.org/records/8279545) - Mesoscale Hydrologic Model.       
* [ORCHIDEE](https://forge.ipsl.jussieu.fr/orchidee/browser/branches/ORCHIDEE-MICT/tags/ORCHIDEE_MICT_8.4.1) - Organising Carbon and Hydrology In Dynamic Ecosystems. 
* [PCR-GLOBWB](https://github.com/UU-Hydro/PCR-GLOBWB_model) - PCRaster Global Water Balance model.
* [SACRAMENTO](https://github.com/NOAA-OWP/sac-sma) - Sacramento Soil Moisture Accounting Model.
* [SWAT](https://swatplus.gitbook.io/docs/source-code) - Soil and Water Assessment Tool.  
* [VIC](https://github.com/UW-Hydro/VIC) - Variable Infiltration Capacity model.

## Replication

We provide all the functions needed to replicate our workflow in the "functions" folder. 

#### Generated data

The "datasets" folder contains the data used in this study. 

* `call_metrics` folder: It contains the function calls per model and language implementation. These
datasets are used to build the call graphs and identify the risky paths in each model.

* `descriptive_statistics` folder: It includes descriptive metrics for each model and
language implementation, such as lines of code, files and numbers of functions.

* `results_per_function` folder: It includes several metrics at the function level.

* `cyclomatic_complexity_functions.csv`. It provides the cyclomatic complexity of
each function in each model and language implementation.

### Functions

The "functions" folder contains all the custom functions coded for the analysis.
They are all sourced from the `.R`, `.pdf` and `.Rmd` files and therefore the 
user of the code does not need to source them separately.

### Code

We offer the code in `.R`, `.pdf` and `.Rmd`. There are three main analyses:

* `code_synthetic_example`: workflow to create the synthetic call graph used to illustrate our
method and check its internal consistency.

* `code_uncertainty_analysis`: workflow to show the potential of implementing uncertainty
and sensitivity analysis to the method.

* `code_hydrological_models`: application of the method to the fourteen global land use
and hydrological models listed above.

Our entire workflow can be run and the results replicated from either of these files. 
The user must run the code from the same folder where the files in the generated data 
section are stored for a successful compilation.


