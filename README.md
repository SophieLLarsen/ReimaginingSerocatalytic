# ReimaginingSerocatalytic
Code and data for: "Reimagining the Serocatalytic Model for Infectious Diseases"

All codes and data necessary to produce simulation data, aggregate and analyse this data, and produce the figures for this paper are included.

## Requirements
This code can be run with R versions 4.2.1-4.4.1. Required libraries are listed/loaded at the top of each code file. 

## Running the code
Code to fit models, run model simulations on best-fit parameters, and run individual stochastic trajectories can be found in "Hypotheses.Rmd". The anticipated run time is ~3-4 hours. Folders starting with "Profile" contain R code and data to run parameter profiling for a given sHCoV and model, on a high-performance computer. Each job runs via a Bash script that parallelizes profiling for all parameters. Please note that this can take up to several hours per parameter. "Main.Rmd" contains code to generate main figures and tables. "Supplement.Rmd" contains supplementary figure and table code. 

## Serological samples
Between January and March 2020, a cross-sectional study was undertaken [1] at Guangzhou and Red Cross, Hong Kong on seroprevalence of four HCoVs in volunteers from Guangdong Women and Children Hospital and The Chinese University of Hong Kong. Plasma samples were collected from 1886 pediatric patients under 18 years old without signs of influenza-like illness as well as 528 volunteers whose age ranging from 16-67 years old. All peripheral blood samples were centrifuged at 3000 x g for 10 minutes at room temperature for plasma collection and kept at -80°C until used. All study procedures were performed after informed consent. The study was approved by the Human Research Ethics Committee at the Guangdong Women and Children Hospital (Approval number: 202101231) and The Chinese University of Hong Kong (IRB: 2020.229).

## Other data
We used serological samples from [2] to parameterize the seroreversion rates in our study. 


## References
[1] Y. Luo, H. Lv, S. Zhao, Y. Sun, C. Liu, C. Chen, W. Liang, K. Kwok, Q. W. Teo, R. T. So, Y. Lin, Y. Deng, B. Li, Z. Dai, J. Zhu, D. Zhang, J. Fernando, N. C. Wu, H. M. Tun, and X. Mu. Age related seroprevalence trajectories of seasonal coronaviruses in children including neonates in guangzhou, china. International Journal of Infectious Diseases, 127:26–32, 2023.390

[2] A.W.D. Edridge, J. Kaczorowska, A.C.R. Hoste, M. Bakker, M. Klein, K. Loens, M.F. Jebbink, A. Matser, C.M. Kinsella, P. Rueda, et al. Seasonal coronavirus protective immunity is short-lasting. Nature medicine, 26(11):1691–1693, 2020.395
