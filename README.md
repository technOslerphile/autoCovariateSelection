# autoCovariateSelection
**R Package to Implement Automated Covariate Selection for Two Exposure Cohorts Using High-Dimensional Propensity Score Algorithm**

[**Click here to go to the online R Package Documentation**](https://technoslerphile.github.io/autoCovariateSelection/index.html)

[**Now available in CRAN**](https://cran.r-project.org/web/packages/autoCovariateSelection/index.html)
As of February 2026, this package has been downloaded more than 20,000 times from CRAN.

**About this package**

*autoCovariateSelection* is an R package that contains functions to implement automated covariate selection using methods described in the high-dimensional propensity score (HDPS) algorithm put forward by Schneeweiss et.al. 

**Article with details on HDPS implementation**

Karim, M. E. (2025). High-Dimensional Propensity Score and Its Machine Learning Extensions in Residual Confounding Control. The American Statistician, 79(1), 72–90. https://doi.org/10.1080/00031305.2024.2368794

**Detailed tutorial on HDPS and its implementation using this R package**

Karim, M. E. (2023a), “High-Dimensional Propensity Score and its Machine Learning Extensions in Residual Confounding Control in Pharmacoepidemiologic Studies,” Zenodo, DOI: 10.5281/zenodo.7877767.https://ehsanx.github.io/hdPSw/

**Reference paper on HDPS**

Schneeweiss S, Rassen JA, Glynn RJ, Avorn J, Mogun H, Brookhart MA. High-dimensional propensity score adjustment in studies of treatment effects using health care claims data. Epidemiology. 2009 Jul;20(4):512-22. doi: 10.1097/EDE.0b013e3181a663cc. Erratum in: Epidemiology. 2018 Nov;29(6):e63-e64.PMID: 19487948; PMCID: PMC3077219. [link](https://doi.org/10.1097/ede.0b013e3181a663cc)

**Notes on this package**

Covariate adjustment in real-world-observational-data (RWD) is important for estimating adjusted outcomes and this can be done by using methods such as, but not limited to, propensity score matching, propensity score weighting and regression analysis. While these methods strive to statistically adjust for confounding, the major challenge is in selecting the potential covariates that can bias the outcomes comparison estimates in observational RWD (Real-World-Data). This is where the utility of automated covariate selection comes in. 

The functions in this package help to implement the three steps of automated covariate selection. These three functions, in order of the steps required to execute automated covariate selection are:
1. get_candidate_covariates() 
2. get_recurrence_covariates()
3. get_prioritised_covariates(). 

In addition to these functions, a sample real-world-data from publicly available de-identified medical claims data is also available for running examples and also for further exploration. 

**How to install the package**
```
install.packages("autoCovariateSelection")
```
**Example using inbuilt data within the package**
```
library("autoCovariateSelection")
library("dplyr")
data(rwd)
basetable <- rwd %>% select(person_id, treatment, outcome_date) %>% distinct()
patientIds <- basetable$person_id #this can be used as patientIdVector argument for the functions in this package
step1 <- get_candidate_covariates(df = rwd,  domainVarname = "domain", eventCodeVarname = "event_code" ,
                                  patientIdVarname = "person_id", patientIdVector = patientIds,n = 100, min_num_patients = 10)
out1 <- step1$covars_data #this will be input to get_recurrence_covariates() function
step2 <- get_recurrence_covariates(df = out1, patientIdVarname = "person_id", eventCodeVarname = "event_code", patientIdVector = patientIds)
out2 <- step2$recurrence_data #this will be input to get_prioritised_covariates() function
step3 <- get_prioritised_covariates(df = out2, patientIdVarname = "person_id", exposureVector = basetable$treatment,
                                    outcomeVector = ifelse(is.na(basetable$outcome_date), 0,1),patientIdVector = patientIds, k = 10)
out3 <- step3$autoselected_covariate_df #dataframe containing the auto selected covariates for all the subjects in the table  
```
--------------------------------------------------------------------------------------------------------
I gladly welcome all suggestions for improvements and collaborations

[*Please report any bugs or issues here*](https://github.com/technOslerphile/autoCovariateSelection/issues)

**Please cite this package if you use it for your research publications**

  Robert D (2020). autoCovariateSelection: Automated Covariate Selection Using HDPS Algorithm. R package version 1.0.0,
  <https://CRAN.R-project.org/package=autoCovariateSelection>.

