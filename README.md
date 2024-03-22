# IFNG_ICB
Datasets and analysis code for study examining the impact of IFNG signaling alterations on immune checkpoint blockade response.
Title: "Alterations in IFN-γ signaling genes sensitize tumors to immune checkpoint blockade"

This repository is structured into 3 folders:  
1) all-vars: all identified variants are included.  
2) pathogenic-vars-ESM1b: only variants predicted to be pathogenic by ESM1b are counted (https://www.biorxiv.org/content/10.1101/2022.08.25.505311v1.full). Samples with only predicted benign variants are otherwise excluded.  
3) pathogenic-vars-EVE: only variants predicted to be pathogenic by EVE are counted (https://evemodel.org/). Samples with only predicted benign or "uncertain" variants are otherwise excluded.  

Each folder contains:  
1) Patient-level dataset after preprocessing for that particular variant filtering strategy.  
2) Summary of response and alteration frequencies in each dataset, with annotation of which datasets were merged or dropped.  
3) Prespecified sample order used for generating forest plots.  
4) Table of regression results for each study, detailing association of alterations in IFN-γ signaling on ICB response after adjustment for tumor mutation burden (TMB).  
5) Code that will fully reproduce the analysis in the manuscript.  
  
In the patient-level dataset files, the columns "AnyAlt", "PredPathogenicAlt_ESM1b", and "PredPathogenicAlt_EVE" are binary variables that indicate whether a given sample has an alteration depending on the specific criteria used.  

Please contact Ryan Chow with any questions: chowr@pennmedicine.upenn.edu
