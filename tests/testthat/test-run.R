context("get_candidate_covariates")

test_that("testing rwd data", {
  data(rwd)
  expect_equal(nrow(rwd), 69333)
  expect_equal(ncol(rwd), 9)
  expect_equal(length(unique(rwd$person_id)), 1000)
})



test_that("testing end-to-end run of automated covariate selection...", {
  data(rwd)
  basetable <- rwd %>% select(person_id, treatment, outcome_date) %>% distinct()
  patientIds <- basetable$person_id #this can be used as patientIdVector argument for the functions in this package
  step1 <- get_candidate_covariates(df = rwd,  domainVarname = "domain", eventCodeVarname = "event_code" ,
  patientIdVarname = "person_id", patientIdVector = patientIds,n = 100, min_num_patients = 10)
  out1 <- step1$covars_data #this will be input to get_recurrence_covariates() function
  step2 <- get_recurrence_covariates(df = out1, patientIdVarname = "person_id", eventCodeVarname = "event_code", patientIdVector = patientIds)
  out2 <- step2$recurrence_data #this will be input to get_prioritised_covariates() function
  out3 <- get_prioritised_covariates(df = out2, patientIdVarname = "person_id", exposureVector = basetable$treatment,
                                     outcomeVector = ifelse(is.na(basetable$outcome_date), 0,1),patientIdVector = patientIds, k = 10)
  expect_error(get_prioritised_covariates(df = out2, patientIdVarname = "person_id", exposureVector = basetable$treatment,
                                     outcomeVector = ifelse(is.na(basetable$outcome_date), 0,1),patientIdVector = patientIds, k = 0))
  out5 <- get_prioritised_covariates(df = out2, patientIdVarname = "person_id", exposureVector = basetable$treatment,
                                     outcomeVector = ifelse(is.na(basetable$outcome_date), 0,1),patientIdVector = patientIds, k = 1000)
})
