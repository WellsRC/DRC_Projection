# The projected burden of COVID-19 in the Democratic Republic of Congo
Chad R. Wells, Jason K Stearns, Pascal Lutumba, Alison P Galvani

Copyright 2020, Chad Wells et al. All rights reserved. Released under the GGNU GENERAL PUBLIC LICENSE v3.

The MATLAB code provided here will run analysis for the modelling portion of the manuscript.

## System of ordinary differential equations
SODENH - Contains the system of ordinary differential equations
## Demographics
DemoDRC - Compresses the contact matrix, the population, and case fatality ratios for the specified age groups
## Calibration of transmission rate
CalcR0NH - Contains the next generation matrix used in calibrating the transmission rate (beta) for a given R_0. The calibration of beta does not consider that severe/critically ill cases enter their home 24 hrs after symptom onset.
## Parameters
ParameterOutput-Returns the paramters need to run the simulations
## Calculation of average R0 and 95% confidence interval
MIDAS_R0 - Calculates the average R0 and the 95% CI and outputs the numerical results for these values of R0
## Figure 1
RunSim- The script runs the simulations  for the results for Figure 1 and generates the figure
