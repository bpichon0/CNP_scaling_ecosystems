# Code for the CNP lake ecosystem model

> [!NOTE]
> This folder contains all codes to reproduce the analyses and the figures of the C:N:P lake/stream paper.
> Contact: Benoît Pichon, **benoit.pichon0@gmail.com**


This repository contains the code used to perform the analyses (`CNP_lake_main.R`) and plot the figures (`Make_figs.R`).
All the code was made on R (*v4.4.1*) and in julia (*v1.7.3*).

<p align="center">
    <img src="https://github.com/bpichon0/CNP_scaling_ecosystems/blob/master/Figures/Fig_model_lake.jpg" width="800">
</p>

## `Installing R & Julia dependencies`

To install the packages needed for the analyses and create the folder architecture used to save the data sets, please load the file `CNP_lake_functions.R` using: 

```R
source(CNP_lake_functions.R)
```

> [!IMPORTANT]  
> By default, the script assumes that Julia is on the computer, and runs all analyses with Julia.
> If julia is not installed, all analyses can still be ran with R with only few changes.
> First, comment the folling lines in `CNP_lake_functions.R`: L21 to L23 and L309 to L393.
> Then replace L4 of `CNP_lake_main.R` with: *Run_with_julia=F* so that all simulations are made with R.

## Organisation of the script

`CNP_lake_main.R` is the main script of analysis and is devided into four different steps.
 The first two correspond to the analysis of the system with only two functional groups (Step 1 = fixers and decomposers) while (Step 2 = non-fixers and decomposers). Then Step 3 logically corresponds to the full model system. Last, Step 4 corresponds to the same analysis of Step 3 but with a different way of computing indirect effects.

`Make_figs.R` produces the figures from the tables created by the `CNP_lake_main.R` script.

All data and simulations outputs are given so the figures can directly be reproduced using the data in `./data/Simulations`




