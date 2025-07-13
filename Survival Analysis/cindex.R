## -----------------------------------------------------------------------------------------------
library(data.table)
library(mlr3)
library(mlr3learners)
library(mlr3measures)
library(ggplot2)
library(pec)
library(survival)
library(dplyr)
library(randomForestSRC)
library(survAUC) # concordance index



## -----------------------------------------------------------------------------------------------

sps <- fread("../data/sps.csv")

# Drop unwanted columns
drop_cols <- c(
  # we want to drop these
  "dbsource", "formulary_drug_cd_list", "admission_location", "admittime", "dischtime",
  "edregtime", "edouttime", "diagnosis", "latest_gcs_time", "comorbidities", "proc.icd9_list",
  
  # dropping but not sure how to handle yet
  "gcs_total", "gcs_verbal", "gcs_motor", "gcs_eye", "language", "drug.Miscellaneous", "diag.Missing"
)

# dropping those columns
sps[, (drop_cols) := NULL]

# drop patients with missing age (until we figure out what imputation we're doing)
sps <- sps[!is.na(age.at.admit)]





## -----------------------------------------------------------------------------------------------
# ELIANA'S SPS MERGE CODE
sps <- sps %>% mutate(ethnicity = case_when(
    str_detect(ethnicity, regex("ASIAN", ignore_case = TRUE)) ~ "Asian",
    str_detect(ethnicity, regex("WHITE|MIDDLE", ignore_case = TRUE)) ~ "White",
    str_detect(ethnicity, regex("BLACK", ignore_case = TRUE)) ~ "Black or African American",
    str_detect(ethnicity, regex("OTHER|PATIENT|UNABLE|UNKNOWN", ignore_case = TRUE)) ~ "Unknown",
    str_detect(ethnicity, regex("AMERICAN", ignore_case = TRUE)) ~ "American Indian or Alaska Native",
    str_detect(ethnicity, regex("HISPANIC", ignore_case = TRUE)) ~ "Hispanic",
    str_detect(ethnicity, regex("MULTI", ignore_case = TRUE)) ~ "More than one race",
    TRUE ~ ethnicity
  ))


## -----------------------------------------------------------------------------------------------
# categorical variables
cat_vars <- c("gender", "admission_type", "insurance", "religion", "marital_status", "ethnicity", "discharge_location", "intervention.group")


# one-hot encode the categorical variables
mm <- model.matrix(~ . - 1, data = sps[, ..cat_vars])

# combining back with the rest of dataset
sps <- cbind(sps[, !cat_vars, with = FALSE], as.data.table(mm))


## -----------------------------------------------------------------------------------------------
# filter out zero survival
sps <- sps[survival_days >= 0]

# Binary target
sps[, survived_90 := survival_days >= 90]

# Fill NAs
cols_to_fill <- setdiff(names(sps), c("age.at.admit", "gcs_total"))
sps[, (cols_to_fill) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = cols_to_fill]

# Remove ID columns
X <- sps[, !c("subject_id", "survival_days", "survived_90", "event"), with = FALSE]
y <- sps[, .(time = survival_days, status = as.integer(event))]

# Combine for mlr3 task
data <- cbind(X, y)


## -----------------------------------------------------------------------------------------------

sps <- sps[, !duplicated(names(sps)), with = FALSE]

# Prepare survival object
sps <- sps %>% dplyr::filter(!is.na(event), !is.na(survival_days))

# Train/test split
set.seed(7)
train_idx <- sample(seq_len(nrow(sps)), size = 0.5 * nrow(sps))
train_data <- sps[train_idx, ]
test_data <- sps[-train_idx, ]

# fit random survival forest
rsf_model <- rfsrc(Surv(survival_days, event) ~ ., data = train_data %>% select(-subject_id, -survived_90))

# fit gradient boosting model (Cox PH boosting)
library(gbm)
gbm_model <- gbm(Surv(survival_days, event) ~ ., 
                 data = train_data %>% select(-subject_id, -survived_90),
                 distribution = "coxph",
                 n.trees = 100)

rsf_pred <- predict(rsf_model, newdata = test_data)

gbm_pred <- predict(gbm_model, newdata = test_data %>% select(-subject_id, -survived_90), n.trees = 100)




## -----------------------------------------------------------------------------------------------
y_test <- Surv(test_data$survival_days, test_data$event)

# risk score supposedly can be approximated as 1 - predicted survival probability at a fixed time (e.g., median), or use predicted mortality (1 - survival) at last time point

# risk score: 1 - survival probability at max time point predicted
risk_scores <- 1 - rsf_pred$survival[, ncol(rsf_pred$survival)]

rsf_pred <- predict(rsf_model, newdata = test_data)
rsf_cindex <- randomForestSRC::get.cindex(time = test_data$survival_days,
                               censoring = test_data$event,
                               predicted = -rsf_pred$predicted)
print(paste("RSF C-index:", round(rsf_cindex, 3)))




## -----------------------------------------------------------------------------------------------
gbm_concordance <- concordance(y_test ~ gbm_pred)
if (gbm_concordance$concordance < 0.5) {
  gbm_concordance$concordance <- 1 - gbm_concordance$concordance
}
print(paste("GBM C-index:", round(gbm_concordance$concordance, 3)))


## -----------------------------------------------------------------------------------------------
time_points <- c(90, 360, 720, 1460)

get_surv_probs <- function(surv_row, time_interest, time_points) {
  sapply(time_points, function(t) {
    idx <- max(which(time_interest <= t))
    surv_row[idx]
  })
}

# printing survival probabilities for first 5 test patients
for (i in 1:5) {
  surv_probs <- get_surv_probs(rsf_pred$survival[i, ], rsf_pred$time.interest, time_points)
  cat(sprintf("Patient %d survival probabilities:\n", i))
  for (j in seq_along(time_points)) {
    cat(sprintf("  P(survival > %d days) = %.3f\n", time_points[j], surv_probs[j]))
  }
}


## -----------------------------------------------------------------------------------------------
gbm_cindex <- Hmisc::rcorr.cens(
  x = -gbm_pred,  # Note the negative sign to ensure proper risk direction
  S = y_test
)["C Index"]
print(paste("GBM C-index (Hmisc):", round(gbm_cindex, 3)))


## -----------------------------------------------------------------------------------------------
rcorr.cens # CORRELATION COEFFICIENT FOR CENSORED DATA


## -----------------------------------------------------------------------------------------------


