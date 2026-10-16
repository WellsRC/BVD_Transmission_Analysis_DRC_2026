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
Returns historical epidemiological parameters for BVD Ebola
### Calculate_Prior_Weights
Calculates the prior weights in the ABC fitting script BVD_ABC_Fitting
### Sample_Approximation_Detection_Chain
A function that approximates the computation of the equation for the detection of cases in a chain
### Sample_BVD_Parameters
Samples the parameters used in the ABC fitting script BVD_ABC_Fitting
### Simulate_Branching_Epidemic
Simulates an epidemic using a branching process
### SIR_ODE
The system of ODEs used to generate the prevalence that serves as a proxy for the seeding rate in the generation of Figure 2

## Script Summary
### BVD_ABC_Fitting
A script used to run the ABC fitting algorithm to the BVD Ebola outbreak in the DRC
### BVD_Epidemic_Analysis_Alternate_mle
Applies the HINT framework under the maximum likelihood parameter estimates for the 2026 BVD outbreak in the DRC to compute the likelihood of the origin date and number of index cases 
### BVD_Epidemic_Analysis
Applies the HINT framework under historical parameter estimates to compute the likelihood of the origin date and number of index cases 
### Epidemic_DRC_Results_Supplementary_Figure_5
Generates the output and figures for the analysis conducted in BVD_Epidemic_Analysis (Supplementary Figure 5)
### Epidemic_DRC_Results_Supplementary_Figure_6_Alternate_mle
Generates the output and figures for the analysis conducted in BVD_Epidemic_Analysis_Alternate_mle (Supplementary Figure 6)
### Figure_1
Generates the results for the fitting of the model to the reported deaths and cases in the epidemic and the origin date estimated from the ABC fitting(Figure 1)
### Figure_2
Generates the results for attribution of cases with unknown exposure to local transmission or seeding. (Figure 2)
### Figure_3
Generates the results for the distribution of transmission-chain detection by generation (Figure 3)
### Results_False_Extinction
Outputs the results of the main text for false extinction
### Supplementary_Figure_1
Outputs the results of the main text for the summary statistics of the parameters estimated during the ABC fitting (Supplementary Figure 1)
### Supplementary_Figure_2_and_3 
Generates the figures computing the odds ratio of localized transmission relative to a seeding event across basic reproduction numbers for a homogeneous seeding rate (Supplementary Figures 2 and 3)
### Supplementary_Figure_4
Generates the figures computing the PRCC values (Supplementary Figure 4)

## File Summary
### BVD_Model_Fit
The output from the ABC fitting script BVD_ABC_Fitting
### Epidemic_Analysis_Likelihood_Alternate_mle
The output from the script BVD_Epidemic_Analysis_Alternate_mle
### Epidemic_Analysis_Likelihood
The output from the script BVD_Epidemic_Analysis
### Extinction_Detection_Analysis
The output save from the initial part of Figure_3 to reduce the time in regenerating the Figure
### Output_Fig_S4
The output saved from the initial part of Supplementary_Figure_4 to reduce the time in regenerating the Figure
