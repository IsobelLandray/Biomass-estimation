# Biomass-estimation
Update on estimation of global human biomass

`Height_NCDRiSC.Rmd` cleans the height data from NCD-RisC and outputs "meanheight_2022.rds". 

`DATASET_cleaning.Rmd` then cleans the BMI data and height data from NCD-RisC and merges with UN population data to output "df_bmipopheight_2022_18plus.csv" which is the main dataset used for the analysis

The R file `halfnormals_pivotprior_implementation_PLUSNORMAL.R` implements the Stan files for the split-normal and normal distribution Bayesian multinomial models `halfnormals_scaled_pivotprior.stan` and `normal_forcompar.stan`, respectively. The datase "twohalfnorms_results_newmid_rules_modeshift_2022_31Oct_pivotprior_pluslimit_plusnormal_updatedpriors.RDS" is the dataset with the posterior estimates by age-sex-country group from the Bayesian multinomial model with split-normal distribution and with normal distribution for BMI (further details on these models, with different age groups are detailed in the code in the repository: https://github.com/IsobelLandray/Modelling_BMI/tree/main).

`WON_analysiswith2022data_16April_excessbiomassversion.qmd` is where all the analysis for the manuscript is conducted: including figures, tables, results in text. "Supplementary_Table1_resultsbycountryagesexgroup_17Jun.rds" is the Supplementary Table 1 from the paper which gives BMI distribution parameter estimates and biomass calculations under the three scenarios (current, lowest, highest) for each age-sex-country group.
