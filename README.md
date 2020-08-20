This repository contains the code and data to reproduce the figures and major analyses in 

Grohn J, Schüffelgen U, Neubert F-X, Verhagen L, Sallet J, Kolling N, Rushworth MFS. Multiple systems in macaques for tracking prediction errors and other types of surprise. PLOS Biology. 2020.

The `main.m` script generates the figures by loading in the `data` folder and executing the scripts in the `scripts` folder. The `timeseries` folder contains extracted BOLD timecourses that are also loaded in and plotted by the `main.m` script. Finally, the `maps` folder contains non-cluster-corrected whole-brain maps of the z-statistics described in the paper. They can be viewed with [fsleyes](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLeyes) or similar software. 
