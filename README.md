# missing-values-diagnostic-studies
Code for the simulation study and case study 

This project belongs to the manuscript entitled "Comparison of methods to handle missing values in a continuous index test in a diagnostic accuracy study – a simulation study" by K Stahlmann (1), B Kellerhuis (2), JB Reitsma (2), N Dendukuri (3), A Zapf (1). 
1 Institute of Medical Biometry and Epidemiology, University Medical Center Hamburg Eppendorf, Hamburg, Germany
2 Julius Center for Health Sciences and Primary Care, University Medical Center Utrecht, Utrecht University, Utrecht, The Netherlands
3 Department of Medicine, McGill University, Montreal, Canada

Author of the code: Katharina Stahlmann (k.stahlmann@uke.de)
This project includes the code for the simulation study (folder Analyse) and case study (folder Case study) including intermittent results. 

The folder Analyse contains all scripts and intermittent results used within the simulation study. The RMarkdown files "Simulation_master.Rmd", "Simulation_supplement1.Rmd" and "Simulation_supplement2_(mnar).Rmd" are used for analysing the simulated AUC data and produce the results of the main analysis, supplemental additional methods and supplemental MNAR scenarios (corresponding to the additional files Additional_file_4 (including the figures in the mansucript), Additional_file_1 and Additional_file_2.html). The files "Simulation_performance3.R" and "mnar_test.RData" are used within these files. As several methods are very time consumig, we had to execute our simulation on a high performance cluster divided in parts and separately for (nearly) each method. This resulted in many files. All files for the simulation can be found under Analyse/hpc and the results of the simulation under Analyse/hpc/Simulation_data. The latter are used by the Rmarkdown files to produce the results presented in the manuscript. The simulation was conducted with R Version 4.1.0 on a high performance cluster running on CentOS 7 Linux (operated by the University Hamburg). The simulated data were analysed on a local maching with R Version 4.4.1. 

The folder Case study contains all scripts and data used within the case study. The RMarkdown file "Case_study_analysis.Rmd" analyses the data (nnhs2.dta) under use of the R scripts Calculations.R and functions_AUC_ROC.R. It produces the results of Additional_file_3.html and the figures in the manuscript. Internittent results of the case study analysis can be accessed under /Case study/Ergebnisse. The case study was analysed on a local maching with R Version 4.4.1. The data nnhs2.dta were obtained from the Diagnostic and Biomarkers Statistical (DABS) Center (Fred Hutch Cancer Center. https://research.fredhutch.org/diagnostic-biomarkers-center/en/datasets.html. Accessed 21 Dec 2023)

If someone is interested in the functions of the methods which are compared in this study, these can be found in Analyse/hpc/functions_AUC.R (the methods produce only the AUC and its confidence interval) and Case study/functions_AUC_ROC.R (AUC+CI and the ROC are produced). These functions are based on the code of the original authors of the respective method. Their references are given within our mansucript, the Additional files and within the functions_AUC.R and functions_AUC_ROC.R files. 

The sessionInfo() for each RMarkdown file can be found in the respective Additional file (Additional_file_1 to _4). 
