# autoCovariateSelection
**R Package to Implement Automated Covariate Selection for Two Exposure Cohorts Using High-Dimensional Propensity Score Algorithm**

[**Click here to go to the online R Package Documentation**](https://technoslerphile.github.io/autoCovariateSelection/index.html)

[**Now available in CRAN**](https://cran.r-project.org/web/packages/autoCovariateSelection/index.html)

*autoCovariateSelection* is an R package that contains functions to implement automated covariate selection using methods described in the high-dimensional propensity score (HDPS) algorithm put forward by Schneeweiss et.al. 

**Original publication by [Schneeweiss et.al. (2009)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3077219/)**

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
