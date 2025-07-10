# Code for "MIRAGE: a Bayesian statistical method for gene-level rare variant analysis incorporating functional annotations"

This repository contains code and data to reproduce the results in our paper:
MIRAGE: a Bayesian statistical method for gene-level rare variant analysis
incorporating functional annotations. 

## Requirements
- R version 4.2+ 
- Packages required `SKAT`, `ACAT`, `AssotesteR`

## Running the analysis

+ `Simulations.R` in folder code is the main code to generate simulated data, run EM algorithm and calculate p values by other methods.
+ `code/Figure S2.R` is the R code to draw Figures in S2, with randm number of variants in a genes when simulating the data.
+ `code/Figure 2 and S2.R` is the code to draw Figure 2 and S2 with same number of variants in a gene when simulating the data.
+ `code/Figure 2 and S4.R` is the code to draw sensitivity figures in Figure 2 and S4. 

## R package 

This is the R package [MIRAGE](https://xinhe-lab.github.io/mirage/). 


## Citation
