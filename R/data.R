#' Sample Data for autoCovariateSelection
#'
#' This is data contains Medicare claims data of a small sample of 1000 patients from the publicly available CMS Medicare De-SynPUF data. It contains
#' all data from three domains - diagnosis, procedures and medications. The diagnosis codes are ICD9 codes, procedures are CPT4/HCPCS codes and medications
#' are NDC codes.
#'
#' @format A data frame with 69333 rows and 9 variables:
#' \describe{
#'   \item{person_id}{patient_identifier}
#'   \item{index_date}{Date of first exposure. For one patient, there will only be one index_date}
#'   \item{event_date}{Date at which event_code occurred for the patient}
#'   \item{event_code}{The medical coding of the event. These are ICD9, CPT4, HCPCS or NDC codes depending on the \code{domain}}
#'   \item{event_concept_id}{Another identifier for the \code{event_code}. This is irrelevant for this package and you can ignore it}
#'   \item{domain}{The domain to which the \code{event_code} belongs to. The three unique values are dx (for diagnosis), px (for procedure)
#'   and rx (for medication)}
#'   \item{treatment}{Binary indicator treatment allocation based on exposure. 1 indicates primary cohort and 0 for control/comparator cohort}
#'   \item{outcome_date}{Date in which the outcome occurred. \code{NA} indicates no outcome occurred. In this sample data, the outcome is death}
#'   \item{last_enrollment_date}{Last enrolled date of the patient. This field is irrelevant for this package and you can ignore it}
#'   ...
#' }
"rwd"
