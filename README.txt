## Paper (under review) 

E. Cookson and R.L. Detwiler (2022). Global patterns and temporal trends of perfluoroalkyl substances (PFAS) in Wastewater.  Manuscript submitted for publication to Water Research


## Overview

We present a meta-analysis of PFAS in wastewater reported in 44 peer-reviewed publications that include 440 influent and 508 effluent samples, collected from 21 countries, for which some or all of five PFCAs (PFHxA, PFHpA, PFOA, PFNA, PFDA) and three PFSA (PFBS, PFHxS, PFOS) were measured. 

We use a linear mixed effect regression model of PFAS concentrations measured in wastewater effluent (collected from 2004 to 2020) to determine temporal trend of global wastewater effluent concentrations of each of the PFAS and the corresponding mean concentration for each country. 

In addition, we compare PFAS concentrations between (i) influent and effluent samples, (ii) samples with different source types, and (iii) liquid and solid suspended particulate matter samples collected from wastewater treatment plants. 


Note: Because we are still waiting approval to publish all data, results may vary slightly from the paper in review. 

### data

used for temporal regression and analysis of GDP per capita 

------ SI_T1.xlsx ------	
Includes the entire meta-dataset of PFAS measured in wastewater (excluding data from Coggan et al., 2019). 

------ GDP_per_capita.xlsx ------	
Includes GDP per capita data for all countries analyzed in the temporal regression.  
(The World Bank, 2020. GDP per capita., https://data.worldbank.org/indicator/NY.GDP.PCAP.CD. Accessed: 2021-09-27)


used for correlating PFAS or comparing between samples 

------ SI_T2.xlsx ------	
Meta-dataset of PFAS measured in wastewater. Includes only reported mean concentrations or measurements of PFAS that could be connected to a specific sample. Because the temporal regression analysis considered each PFAS independently, we could use all reported data for that analysis. 



### code

------ TemporalRegression.m (run) ------
Performs linear mixed effect regression on log PFAS effluent concentrations and determines mean concentrations for each country.

------ plotTemporalRegression.m (function) ------
Plots observations, linear model, and linear mixed effect model regression results with respect to sample time. 

------ plot_bprime.m (function) ------
Generates figure of temporal regression results (intercepts by country, mean intercept of PFAS, and slopes of PFAS) 

------ GDPpercapRegression.m (function) ------
Performs linear regression on observations shifted to 2019 (by linear mixed effect model slope) on national 2019 GDP per capita and generates figure.


------ inf_vs_eff (run) ------
Plots effluent concentrations vs influent concentrations

------ SourceType_PCA.m (run)------
Principal component analysis of PFAS effluent concentrations, plots with respect to wastewater source type

------ cell_str_2_num.m (function) ------ 
Converts ND values from '<LOD' to 0.5*LOD	 

