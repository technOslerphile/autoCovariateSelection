#' Generate candidate empirical baseline covariates based on prevalence in the baseline period
#' @import dplyr
#' @description
#' \code{get_candidate_covariates} function generates the list of candidate empirical covariates based on their prevalence
#' within each domains (dimensions). This is the first step in the automated covariate selection process. See 'Automated Covariate Selection'
#' section below for more details regarding the overall process.
#'
#' @section Automated Covariate Selection:
#' \strong{The three steps in automated covariate selection are listed below with the functions implementing the methodology}
#' \enumerate{
#' \item Identify candidate empirical covariates: \code{\link[autoCovariateSelection]{get_candidate_covariates}}
#' \item Assess recurrence: \code{\link[autoCovariateSelection]{get_recurrence_covariates}}
#' \item Prioritize covariates: \code{\link[autoCovariateSelection]{get_prioritised_covariates}}
#' }
#'
#' @author Dennis Robert \email{dennis.robert.nm@gmail.com}
#'
#' @param df The input \code{data.frame}. This should contain at least 3 fields containing information on patient identifier,
#' covariate codes and domain names of covariate codes in a long format. Any other fields containing values such as dates,
#' treatment group are optional and will be ignored for this analysis
#' @param domainVarname The variable(field) name which contains the domain of the covariate in the \code{df}.
#' The domains are usually diagnosis, procedures and medications.
#' @param eventCodeVarname The variable name which contains the covariate codes (eg:- CCS, ICD9) in the \code{df}
#' @param patientIdVarname The variable name which contains the patient identifier in the \code{df}
#' @param patientIdVector The 1-D vector with all the patient identifiers. The length of this vector should be equal to
#' the number of distinct patients in the \code{df}. This vector is not really used in the function analysis per se. This is
#' used only to return the same back as function output because the filtered \code{df} based on \code{covars} will likely not
#' contain all patients in the input \code{df} because there could be patients for whom no records were found for any of the
#' identified \code{covars} and they will thus be not present in the filtered \code{df} which is also an output of this function.
#' The \code{patientIds} vector output will contain all original patients and by returning this vector, it can later be  used in the
#' next steps of automated covariate selection because each step is dependent on previous steps and information on patients who did not
#' have any identified \code{covars} is also important for the next steps. This is why this vector is an input as well as an output, without
#' affecting the analysis of this function.
#' @param n The maximum number of empirical candidate baseline covariates that should be returned within each domain.
#' By default, n is 200
#' @param min_num_patients Minimum number of patients that should be present for each covariate to be selected for selection.
#' To be considered for selection, a covariate should have occurred for a minimum \code{min_num_patients} in the baseline period

#' @return
#' A named list containing three R objects
#' \itemize{
#'   \item \code{covars}   A 1-D vector containing the names of selected baseline covariate names from each domain.
#'   For each domain in the \code{df}, the number of \code{covars} would be equal to or less than \code{n}
#'   \item \code{covars_data}   {The \code{data.frame} that is filtered out of \code{df} with only the selected \code{covars}}. The values of the
#'   \code{eventCodeVarname} field is prefixed with the corresponding \code{domain} name. For example, if the event code is 19900 and the domain
#'   is 'dx', then the the covariate name will be 'dx_19900'.
#'   \item \code{patientIds}   {The list of patient ids present in the original input \code{df}. This is exactly the same as the input \code{patientIdVector}}
#'}
#' @examples
#' library("autoCovariateSelection")
#' data(rwd)
#' head(rwd, 3)
#' #select distinct elements that are unique for each patient - treatment and outcome
#' basetable <- rwd %>% select(person_id, treatment, outcome_date) %>% distinct()
#' head(basetable, 3)
#' patientIds <- basetable$person_id
#' step1 <- get_candidate_covariates(df = rwd,  domainVarname = "domain",
#' eventCodeVarname = "event_code", patientIdVarname = "person_id",
#' patientIdVector = patientIds,n = 100, min_num_patients = 10)
#' out1 <- step1$covars_data #this will be input to get_recurrence_covariates() function

#' @details
#' The theoretical details of the high-dimensional propensity score (HDPS) algorithm is detailed in the publication listed below in the \code{References} section.
#' \code{get_candidate_covariates} is the function implementing what is described in the 'Identify candidate empirical covariates' section
#' of the article.
#'
#' @references
#' Schneeweiss S, Rassen JA, Glynn RJ, Avorn J, Mogun H, Brookhart MA. High-dimensional propensity score adjustment in studies of treatment effects using health care claims data Epidemiology. 2009;20(4):512-522. doi:10.1097/EDE.0b013e3181a663cc
#'
#' @export

get_candidate_covariates <- function(df,
                                     domainVarname,
                                     eventCodeVarname,
                                     patientIdVarname,
                                     patientIdVector,
                                     n = 200,
                                     min_num_patients = 100){
  n_patients <- prevalence <- prevalence_truncated <- NULL #binding variables locally to the function to avoid warnings using CMD checks
  cov <- df[c(patientIdVarname, eventCodeVarname, domainVarname)]
  colnames(cov) <- c("patientIdVarname", "eventCodeVarname", "domainVarname")
  sampleSize <- length(unique(patientIdVector))
  prevalences <- (cov %>% group_by(domainVarname, eventCodeVarname) %>% mutate(n_patients = length(unique(patientIdVarname))) %>%
                    mutate(prevalence = (n_patients/sampleSize)*100) %>%
                    ungroup() %>% select(domainVarname, eventCodeVarname, prevalence, n_patients) %>% distinct() %>%
                    filter(n_patients >=min_num_patients) %>%
                    mutate(prevalence_truncated = ifelse(prevalence > 50, 100-prevalence, prevalence)) %>%
                    arrange(domainVarname, desc(prevalence_truncated)) %>%
                    group_by(domainVarname) %>% mutate(rank = seq(1,length(domainVarname))))
  covs_candidates <- prevalences %>% group_by(domainVarname) %>% filter(rank <= n) %>% select(domainVarname, eventCodeVarname, prevalence) %>% distinct()

  outdf <- (cov %>% inner_join(covs_candidates) %>%
              mutate(eventCodeVarname = paste0(domainVarname, "_", eventCodeVarname)) %>%
              select(patientIdVarname, eventCodeVarname))

  colnames(covs_candidates) <- c(domainVarname, eventCodeVarname, "prevalence")
  colnames(outdf) <- c(patientIdVarname, eventCodeVarname)

  return(list("covars"= covs_candidates, "covars_data" = outdf, "patientIds" = patientIdVector))
}


#' Generate the binary recurrence covariates for the identified candidate empirical covariates
#' @import dplyr
#' @importFrom data.table data.table
#' @importFrom data.table dcast
#' @importFrom purrr reduce
#' @importFrom stats quantile
#' @description
#' \code{get_recurrence_covariates} function assesses the recurrence of each of the identified candidate empirical covariates
#' based on their frequency of occurrence for each patient in the baseline period and generates three binary recurrence covariates
#' for each of the identified candidate empirical covariates. This is the second step in the automated covariate selection process.
#' The first step of identifying empirical candidate covariates is done via \code{\link[autoCovariateSelection]{get_candidate_covariates}} function.
#' See 'Automated Covariate Selection'section below for more details regarding the overall process.
#'
#' @section Automated Covariate Selection:
#' \strong{The three steps in automated covariate selection are listed below with the functions implementing the methodology}
#' \enumerate{
#' \item Identify candidate empirical covariates: \code{\link[autoCovariateSelection]{get_candidate_covariates}}
#' \item Assess recurrence: \code{\link[autoCovariateSelection]{get_recurrence_covariates}}
#' \item Prioritize covariates: \code{\link[autoCovariateSelection]{get_prioritised_covariates}}
#' }
#'
#' @author Dennis Robert \email{dennis.robert.nm@gmail.com}
#'
#' @param df The input \code{data.frame}. Ideally this should be the output \code{covars_data}
#' from \code{\link[autoCovariateSelection]{get_candidate_covariates}}
#' @param patientIdVarname The variable name which contains the patient identifier in the \code{df}
#' @param eventCodeVarname The variable name which contains the covariate codes (eg:- CCS, ICD9) in the \code{df}
#' @param patientIdVector The 1-D vector with all the patient identifiers. This should contain all the patient IDs in the original two
#' cohorts. This vector can simply be the \code{patientIds} output vector of the \code{get_candidate_covariates} function.
#' of the function
#'
#' @return
#' A named list containing two R objects
#' \itemize{
#'   \item \code{recurrence_data}   A \code{data.frame} containing all the binary recurrence covariates for all the patients in wide format.
#'   This means that this \code{data.frame} will have a dimension with number of rows equal to number of distinct patients and number of
#'   columns equal to number of binary recurrence covariates plus 1 (for the patient Id variable). The binary recurrence covariate is prefixed with
#'   a 'rec_' to indicate that the covariate is a 'reccurrence covariate' and suffixed with '_once', '_sporadic' or '_frequent'.
#'   See \code{details} section above for details.
#'   \item \code{patientIds}   {The list of patient ids present in the original input \code{df}. This is exactly the same as the input \code{patientIdVector}}
#'}
#'
#' @examples
#' library("autoCovariateSelection")
#' data(rwd)
#' head(rwd, 3)
#' basetable <- rwd %>% select(person_id, treatment, outcome_date) %>% distinct()
#' head(basetable, 3)
#' patientIds <- basetable$person_id
#' step1 <- get_candidate_covariates(df = rwd,  domainVarname = "domain",
#' eventCodeVarname = "event_code" , patientIdVarname = "person_id",
#' patientIdVector = patientIds,n = 100, min_num_patients = 10)
#' out1 <- step1$covars_data
#' all.equal(patientIds, step1$patientIds) #should return  TRUE
#' step2 <- get_recurrence_covariates(df = out1, patientIdVarname = "person_id",
#' eventCodeVarname = "event_code", patientIdVector = patientIds)
#' out2 <- step2$recurrence_data

#' @details
#' The recurrence covariates are generated based on the frequency (counts) of occurrence of each empirical candidate covariates that got
#' generated by the \code{generate_candidate_covariates} function. This is done by looking at the baseline period of each patients and
#' assessing whether the covariate occurred only once or sporadically or frequently. That is, a maximum of three recurrence covariates
#' for each candidate covariate is created and returned.
#' \itemize{
#' \item \code{once} Indicates whether or not the covariate occurred more than or equal to 1 number of times for the patient
#' \item \code{sporadic} Indicates whether or not the covariate occurred more than or equal to median (median of non-zero occurrences of
#' the candidate covariate) number of times for the patient.
#' \item \code{frequent} Indicates whether or not the covariate occurred more than or equal to upper quartile (75th percentile of non-zero
#' occurrences of the candidate covariate) number of times for the patient
#' }
#' Note that if two or all three covariates are identical for any of the binary recurrence covariates, only the distinct recurrence covariate
#' is returned. For example, if once == sporadic == frequent for the candidate covariate (median and upper quartile both are 1), then only the 'once' recurrence covariate is
#' returned. If once != sporadic == frequent, then 'once' and 'sporadic' is returned. If once == sporadic != frequent, then 'once'
#' and 'frequent' are returned. If none of three recurrence covariates are identical, then all three are returned.
#' The theoretical details of the algorithm implemented is detailed in the publication listed below in the \code{References} section.
#' \code{get_recurrence_covariates} is the function implementing what is described in the 'Assess Recurrence' section
#' of the article.
#'
#' @references
#' Schneeweiss S, Rassen JA, Glynn RJ, Avorn J, Mogun H, Brookhart MA. High-dimensional propensity score adjustment in studies of treatment effects using health care claims data Epidemiology. 2009;20(4):512-522. doi:10.1097/EDE.0b013e3181a663cc
#'
#' @export

get_recurrence_covariates <- function(df, patientIdVarname, eventCodeVarname, patientIdVector){
  #long to wide
  colnames(df) <- c("patientIdVarname", "eventCodeVarname")
  df <- data.table(df)
  baseline_data <- data.table::dcast(df, patientIdVarname ~ eventCodeVarname, value.var = "eventCodeVarname", fun.aggregate = length)

  #add the patients who did not have any events related to candidate covariates in baseline period to the baseline_data
  #these patients will have 0 as value for all the candidate covariates as they did not experience the event
  nullPatients <- setdiff(patientIdVector, baseline_data$patientIdVarname)
  m <- matrix(nrow = length(nullPatients), ncol = ncol(baseline_data))
  dim(m)
  m[,1] <- nullPatients
  m[,2:ncol(m)] <- 0 #assign 0 as value of all the baseline covariates for these patients who did not have any event
  colnames(m) <- colnames(baseline_data)
  baseline_data <- rbind(baseline_data, m)
  recs <- list()
  ids <- baseline_data$patientIdVarname
  baseline_data <- data.frame(baseline_data)
  for (i in 2:ncol(baseline_data)){
    covname <- colnames(baseline_data)[i]
    vals <- as.integer(baseline_data[,covname])
    median <- as.integer(median(vals[vals!=0]))
    uq <- as.integer(quantile(vals[vals!=0], 0.75))

    if (1 == median & median == uq){
      recs[[i-1]] <- data.frame("person_id" = ids,
                                "once" =  as.integer(ifelse(vals >=1, 1, 0)))
      colnames(recs[[i-1]]) <- c("person_id",paste0("rec_",covname, "_once"))

    }else if (1 == median & median != uq){
      recs[[i-1]] <- data.frame("person_id" = ids,
                                "once" =  as.integer(ifelse(vals >=1, 1, 0)),
                                "frequent" = as.integer(ifelse(vals >=uq, 1, 0)))
      colnames(recs[[i-1]]) <- c("person_id",paste0("rec_",covname, "_once"), paste0("rec_",covname, "_frequent"))

    } else if (1 != median & median == uq){
      recs[[i-1]] <- data.frame("person_id" = ids,
                                "once" =  as.integer(ifelse(vals >=1, 1, 0)),
                                "sporadic" = as.integer(ifelse(vals >=median, 1, 0)))
      colnames(recs[[i-1]]) <- c("person_id",paste0("rec_",covname, "_once"), paste0("rec_",covname, "_sporadic"))


    }else {
      recs[[i-1]] <- data.frame("person_id" = ids,
                                "once" =  as.integer(ifelse(vals >=1, 1, 0)),
                                "sporadic" = as.integer(ifelse(vals >=median, 1, 0)),
                                "frequent" = as.integer(ifelse(vals >=uq, 1, 0)))
      colnames(recs[[i-1]]) <- c("person_id",paste0("rec_",covname, "_once"), paste0("rec_",covname, "_sporadic"), paste0("rec_",covname, "_frequent"))

    }
  }
  if(length(unique(df$eventCodeVarname)) != length(recs)) stop("The number of covariates for which reccurrence variables are created does not match the original number of input baseline covariates...")

  #generate a data.frame with all the recurrence covariates
  recurrences <-recs %>% purrr::reduce(inner_join, by="person_id")
  if(nrow(recurrences) != length(patientIdVector)) stop("Number of patients for which recurrence covariates are created does not match the original number of patients")
  recurrences <- recurrences[match(patientIdVector, recurrences$person_id),]
  colnames(recurrences)[1] <- patientIdVarname
  return(list("recurrence_data" = recurrences, "patientIds" = patientIdVector))
}


#' Generate the prioritised covariates from the global list of binary recurrence covariates using multiplicative bias ranking
#' @import dplyr
#'
#' @description
#' \code{get_prioritised_covariates} function assesses the recurrence of each of the identified candidate empirical covariates
#' based on their frequency of occurrence for each patient in the baseline period and generates three binary recurrence covariates
#' for each of the identified candidate empirical covariates. This is the third and final step in the automated covariate selection process.
#' The previous step of assessing recurrence and generating the binary recurrence covariates is done
#' using the \code{\link[autoCovariateSelection]{get_recurrence_covariates}} function.
#' See 'Automated Covariate Selection'section below for more details regarding the overall process.
#'
#' @section Automated Covariate Selection:
#' \strong{The three steps in automated covariate selection are listed below with the functions implementing the methodology}
#' \enumerate{
#' \item Identify candidate empirical covariates: \code{\link[autoCovariateSelection]{get_candidate_covariates}}
#' \item Assess recurrence: \code{\link[autoCovariateSelection]{get_recurrence_covariates}}
#' \item Prioritize covariates: \code{\link[autoCovariateSelection]{get_prioritised_covariates}}
#' }
#'
#' @author Dennis Robert \email{dennis.robert.nm@gmail.com}
#' @param df The input \code{data.frame}. Ideally this should be the output \code{recurrence_data} from the
#' \code{\link[autoCovariateSelection]{get_recurrence_covariates}} function
#' @param patientIdVarname The variable name which contains the patient identifier in the \code{df}
#' @param exposureVector The 1-D exposure (treatment/intervention) vector. The length of this vector should be equal to that of
#' \code{patientIdVector} and \code{outcomeVector}. Also, this should be a binary vector with value of 1 for patients primary cohort 1 and 0 for
#' those in comparator cohort. The order of this vector should resonate the order of patients in \code{outcomeVector} and \code{patientIdVector}
#' @param outcomeVector The 1-D outcome vector indicating whether or not the patient experienced the outcome of interest (value = 1) or not (value =0).
#' The length of this vector should be equal to that of \code{patientIdVector} and \code{exposureVector}. The order of elements in this vector should
#' resonate with the order of patients in \code{exposureVector} and \code{patientIdVector}
#' @param patientIdVector The 1-D vector with all the patient identifiers. This should contain all the patient IDs in the original two
#' cohorts with its length and order equal to and resonating with that of \code{exposureVector} and \code{outcomeVector}
#' @param k The maximum number of prioritised covariates that should be returned by the function. By default, this is 500 as described in the original paper
#'
#' @return
#' A named list containing two R objects
#' \itemize{
#'   \item \code{autoselected_covariate_df} A \code{data.frame} in wide format containing the auto-selected prioritised covariates and their values (1 or 0)
#' for each patients
#'   \item \code{multiplicative_bias}{The absolute log of the multiplicative bias term for each of the auto-selected prioritised covariates}
#'}
#'
#' @examples
#' library("autoCovariateSelection")
#' data(rwd)
#' head(rwd, 3)
#' basetable <- rwd %>% select(person_id, treatment, outcome_date) %>% distinct()
#' head(basetable, 3)
#' patientIds <- basetable$person_id
#' step1 <- get_candidate_covariates(df = rwd,  domainVarname = "domain",
#' eventCodeVarname = "event_code" , patientIdVarname = "person_id",
#' patientIdVector = patientIds,n = 100, min_num_patients = 10)
#' out1 <- step1$covars_data
#' all.equal(patientIds, step1$patientIds) #should be TRUE
#' step2 <- get_recurrence_covariates(df = out1,
#' patientIdVarname = "person_id", eventCodeVarname = "event_code",
#' patientIdVector = patientIds)
#' out2 <- step2$recurrence_data
#' out3 <- get_prioritised_covariates(df = out2,
#' patientIdVarname = "person_id", exposureVector = basetable$treatment,
#' outcomeVector = ifelse(is.na(basetable$outcome_date), 0,1),
#' patientIdVector = patientIds, k = 10)
#'
#' @details
#' To prioritise covariates across data dimensions (domains) should be assessed by their potential for controlling confounding that is not conditional
#' on exposure and other covariates. This means that the association of the covariates with the outcomes (relative risk) should be taken into
#' consideration for quantifying the 'potential' for confounding. Relative risk weighted by the ratio of prevalence of the covariates between the
#' two exposure groups is known as multiplicative bias. The other way to do this would be to use the absolute risk and this would have been the rather
#' straight-forward procedure to quantify the potential for confounding. However, this method would invariably down-weight the association between the
#' covariate and the outcome if the outcome prevalence is small and the exposure prevalence is high which is a common phenomenon seen with comparative
#' effective research using real-world-data by retrospective cohort studies. The multiplicative bias term balances this and generates a quantity for each
#' covariate that is reflective of its confounding potential. By ranking the multiplicative bias, the objective is to choose the top \code{k} number of
#' covariates from this procedure. \code{k}, by default, is 500 as described in the original paper. For further theoretical details of the
#' algorithm please refer to the original article listed below in the \code{References} section. \code{get_recurrence_covariates} is the function
#' implementing what is described in the 'Prioritise Covariates' section of the article.
#'
#' @references
#' Schneeweiss S, Rassen JA, Glynn RJ, Avorn J, Mogun H, Brookhart MA. High-dimensional propensity score adjustment in studies of treatment effects using health care claims data Epidemiology. 2009;20(4):512-522. doi:10.1097/EDE.0b013e3181a663cc
#'
#' @export


get_prioritised_covariates <- function(df, patientIdVarname, exposureVector, outcomeVector, patientIdVector, k = 500){
  if((nrow(df)!=length(exposureVector) ) | (nrow(df) != length(outcomeVector))) stop("Mismatch in the length of either outcomevector or exposurevector with the number of patients in df...")
  basetable <- data.frame("patientIdVector" = patientIdVector, "exposureVector" = exposureVector, "outcomeVector" = outcomeVector)
  if(sum(is.na(basetable))!=0) stop("Missing values in either patient ids/exposurevector/outcomevector...")
  if(k==0) stop("k should atleast 1")


  outcome <- outcomeVector
  treatment <- as.factor(exposureVector)
  covars<- df[,2:ncol(df)]
  covar_prev <- by(covars, treatment, colMeans)

  prev0 <- covar_prev[[levels(treatment)[1]]] #prevalence of covariates in unexposed
  prev1 <- covar_prev[[levels(treatment)[2]]] #prevalence of covariates in exposed

  relativeRisk <- get_relative_risk(df = df, outcomeVec = outcome)
  if (length(relativeRisk) != length(prev0) & length(relativeRisk) != length(prev1)) stop("Mismatch in lengths of relative risk and prevalence...")
  numerator <- (prev1 * (relativeRisk-1) + 1)
  denominator <- (prev0 * (relativeRisk-1)+ 1)
  absLogMB <- abs(log(numerator/denominator)) #multiplicative bias, absolute log
  colNums <- order(absLogMB, decreasing = TRUE, na.last = NA) #rank the covariates based on absolute log of the multiplicative bias term
  attr(colNums, "multiplicative_bias") <- absLogMB[colNums]
  auto_selected_covariates <- colnames(covars)[colNums][1:k]
  auto_selected_covariates <- auto_selected_covariates[!is.na(auto_selected_covariates)] #keep only non-NA values

  autoselected_covariates <- df[c(patientIdVarname,auto_selected_covariates)]
  finalQCcondition = ( (nrow(autoselected_covariates) == length(patientIdVector)) && (ncol(autoselected_covariates)-1 <= k)
                       && ncol(autoselected_covariates) == length(auto_selected_covariates)+1)
  if(!finalQCcondition){stop("Final QC failed...")} else {"Autoselected baseline covariates successfully!"}
  return(list("autoselected_covariate_df" = autoselected_covariates, "multiplicative_bias" = absLogMB))

}



#' Compute relative risk for each of the covariates with respect to outcomes occurred
#' @import dplyr
#'
#' @description
#' \code{get_relative_risk} function is a helper function used within the \code{\link[autoCovariateSelection]{get_prioritised_covariates}} function.
#' This function computes the prevalence in the exposed and that in the unexposed and simply returns the relative risk for all the covariates in the
#' input \code{data.frame}
#'
#' @author Dennis Robert \email{dennis.robert.nm@gmail.com}
#' @param df The input \code{data.frame}. Ideally this should be the output \code{recurrence_data} from the
#' \code{\link[autoCovariateSelection]{get_recurrence_covariates}} function. The first column should be the patient identifier column and
#' all other columns should be binary covariates. The values of these binary columns should be 1 indicating occurrence of covariate and 0 indicating
#' no occurrence of the covariate.
#' @param outcomeVec The 1-D outcome vector indicating whether or not the patient experienced the outcome of interest (value = 1) or not (value =0).
#' The length of this vector should be equal to the number of rows of \code{df}. The order of elements in this vector should resonate with
#' the order of patients in \code{df}
#' @return
#' A 1-D vector containing relative risk of the association between the covariate (confounder) and the outcome. Thus, the length of this vector
#' will be equal to the number of covariates.
#'
#' @export

get_relative_risk <- function(df, outcomeVec){
  if (nrow(df) != length(outcomeVec)) stop ("Mismatch between number of patients and length of outcome vector...")
  rr <- as.numeric()
  for (i in 2:ncol(df)){
    nums <- table(df[,i])
    num0 <- as.numeric(nums["0"])
    num1 <- as.numeric(nums["1"])
    exp_out_table <- table(outcomeVec, df[,i])
    out0 <- as.numeric(exp_out_table[2,1])
    out1 <- as.numeric(exp_out_table[2,2])
    prev0 = out0/num0 #prevalence of outcome in the unexposed
    prev1 = out1/num1 #prevalence of outcome in the exposed
    relative_risk <- prev1/prev0 #the relative risk
    rr[i-1] <- relative_risk
  }
  return(rr)
}
