[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17962642.svg)](https://doi.org/10.5281/zenodo.17962642)
# The topology of software risk in scientific models

[Arnald Puy](https://www.arnaldpuy.com/), Federico Díaz, Ulrike Proske, Olivia Richards,
Seth N. Linga, Samuel Flinders, Carmen Aguiló-Rivera, and Fernando G. Tinetti.

This study proposes a framework to identify risky paths in scientific models; that is, 
sequences of function calls whose potential failure can have a larger cascading effect into
other parts of the software. We illustrate our method with a synthetic example and by applying it
to fourteen global land use and hydrological models, listed in the Models section below.

## Abstract

*Software risk is studied in engineering and computer science but remains largely unexamined in scientific computing. Existing approaches typically assess risk at the level of individual functions or modules and are tied to a particular conception of risk, despite evidence of failures arising along execution paths and risk admitting multiple conceptualizations. Here we present a network-based framework that combines software-quality metrics with call graphs to identify execution paths concentrating systemic risk. The method treats the definition of risk as uncertain by comparing path rankings across a continuum space of risk definitions, thus accommodating the diverse priorities that shape risk assessment in scientific modelling. By illustrating its potential on fourteen models spanning ecology, environmental sciences, hydrology and climate change, we reveal that risk concentrates in a small subset of paths. Our approach provides a transparent basis for guiding testing, refactoring and verification in scientific software.*

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

We offer the code in `.R`, `.pdf` and `.Rmd`. There are four main analyses:

* `code_synthetic_example`: workflow to create the synthetic call graph used to illustrate our
method and check its internal consistency.

* `code_uncertainty_analysis`: workflow to show the potential of implementing uncertainty
and sensitivity analysis to the method.

* `code_hydrological_models`: application of the method to the fourteen global land use, 
hydrological, ecological and climate models listed above.

* `scalability_test`: stress-test of the scalability of our method across call graphs
with different degrees of complexity.

Our entire workflow can be run and the results replicated from either of these files. 
The user must run the code from the same folder where the files in the generated data 
section are stored for a successful compilation.

## Citation

If you use this workflow, please cite:

Puy, A. et al. (2026). Code and Datasets of The Topology of Software Risk in 
Scientific Models. Zenodo. doi:10.5281/zenodo.17962642.

## License

MIT License

Copyright (c) 2026 Arnald Puy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
