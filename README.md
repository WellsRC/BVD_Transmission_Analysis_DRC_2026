# Quantifying unobserved transmission chains and false extinction in emerging epidemics: Application to the 2026 Bundibugyo Ebola outbreak 
Chad R. Wells, Abhishek Pandey, Seyed M. Moghadas, Martial L. Ndeffo-Mbah, Lilia Potter-Schwartz, Lamia Ayaz, Meagan C. Fitzpatrick, Alison P. Galvani

Copyright (C) <2026>, Chad R. Wells et. al. All rights reserved. Released under the GNU General Public License (GPL)

This repository contains codes and data used to Hidden-chain Inference of Transmission (HINT) branching-process framework to quantify false extinction arising from unobserved transmission and distinguish local transmission from external seeding for cases with unknown exposure.

The model code is written in MATLAB and results are saved as MATLAB data files (extension .mat), with plots also being constructed in MATLAB. 

## OS System requirements
The codes developed here are tested on Windows operating system (Windows 11 Pro: 64-bit). However as Matlab is available for most operating systems, codes should run on Mac OSX and Linux as well.

## Installation guide
### MATLAB
Installation instruction for MATLAB can be found at https://www.mathworks.com/help/install/install-products.html. Typical install time for MATLAB on a "normal" desktop is around 30-40 minutes. The current codes were developed and tested on MATLAB R2025b.

## Demo
Figure_1 generates Figure 1 in the main text. All mat files to run this script are avaialble. 

## Instructions for use
To generate the Figures and output of the calculations, select a script from Figures section to run in MATLAB and enter the name in the command line. All mat file are availble to generate figures and conduct the calculations. To run analysis on a different set of parameters, adjust the parameters in the script and enter the name of the script in the command line to run.

## Function Summary
### BVD_Parameters
### BVD_Probability_of_False_Extinction
### Calculate_Prior_Weights
### Compute_Reproduction_Number
### Generation_Time
### NB_Branching_Process
### Parameter_Generation_Time
### Probability_Asymptomatic_Index_Symptomatic_Generation
### Probability_Outbreak_Exticnt_Within_Generation
### Recursive_Summation_Asymptomatic_Chain
### Recursive_Summation_Detection_Chain
### Sample_Approximation_Detection_Chain
### Sample_BVD_Parameters
### Simulate_Branching_Epidemic
### SIR_ODE_Alt
### SIR_ODE

## Script Summary
### BVD_ABC_Fitting
### BVD_Epidemic_Analysis_Alternate_mle
### BVD_Epidemic_Analysis
### BVD_Probability_False_Extinction_Events
### BVD_Probability_of_Extinction
### Epidemic_DRC_Results_Supplementary_Figure_6
### Epidemic_DRC_Results_Supplementary_Figure_7_Alternate_mle
### Figure_1
### Figure_2
### Figure_3
### Results_False_Extinction
### Supplementary_Figure_1
### Supplementary_Figure_2
### Supplementary_Figure_3_and_4
### Supplementary_Figure_5
### Supplementary_Figure_Extinction
### Unknown_Exposure_Spillover_vs_Transmission_Homogeneous_Rate

## File Summary
### BVD_Model_Fit
### Epidemic_Analysis_Likelihood_Alternate_mle
### Epidemic_Analysis_Likelihood
### Ex_vs_Det_Test
### Extinction_Detection_Analysis
### Numerical_Branching_Simulations
### Output_Fig_S5
### Output_Test
