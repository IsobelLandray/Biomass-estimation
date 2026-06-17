# Biomass-estimation
Update on estimation of global human biomass

Height_NCDRiSC.Rmd cleans the height data from NCD-RisC and outputs meanheight_2022.rds. 

twohalfnorms_results_newmid_rules_modeshift_2022_31Oct_pivotprior_pluslimit_plusnormal_updatedpriors.RDS is the dataset with the posterior estimates by age-sex-country group from the Bayesian multinomial model with split-normal distribution for the BMI distribution (see for more details: ).

DATASET_cleaning.Rmd then cleans the BMI data and height data from NCD-RisC and merges with UN population data to output df_bmipopheight_2022_18plus.csv which is the main dataset used for the analysis. 

WON_analysiswith2022data_16April_excessbiomassversion.qmd is where all the analysis for the manuscript is conducted. 

Supplementary_Table1_resultsbycountryagesexgroup_17Jun.rds is the Supplementary Table 1 from the paper which gives BMI distribution parameter estimates and biomass calculations under the three scenarios (current, lowest, highest) for each age-sex-country group.
