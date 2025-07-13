### ROI before imputation 



                            #### HANRUI's CODE START####



rm(list = ls())
library(dplyr)
library(mice)
library(stringr)

df_bf_impute = read.csv("C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/merged_v4.csv")


df_bf_impute <- df_bf_impute %>%
  mutate(ethnicity = case_when(
    str_detect(ethnicity, regex("ASIAN", ignore_case = TRUE)) ~ "Asian",
    str_detect(ethnicity, regex("WHITE|MIDDLE", ignore_case = TRUE)) ~ "White",
    str_detect(ethnicity, regex("BLACK", ignore_case = TRUE)) ~ "Black or African American",
    str_detect(ethnicity, regex("OTHER|PATIENT|UNABLE|UNKNOWN", ignore_case = TRUE)) ~ "Unknown",
    str_detect(ethnicity, regex("AMERICAN", ignore_case = TRUE)) ~ "American Indian or Alaska Native",
    str_detect(ethnicity, regex("HISPANIC", ignore_case = TRUE)) ~ "Hispanic",
    str_detect(ethnicity, regex("MULTI", ignore_case = TRUE)) ~ "More than one race",
    TRUE ~ ethnicity
))





### mice for imputation: (survival days not included)
df_bf_impute = df_bf_impute %>%
  select(-1)
df_bf_impute <- df_bf_impute %>% mutate(across(where(is.character), factor))
meth <- make.method(df_bf_impute)
meth["age.at.admit"] <- "pmm"
meth["gcs_total"]    <- "pmm"
multicat <- c("admission_type", "admission_location", "discharge_location",
              "insurance", "language", "religion",
              "marital_status", "ethnicity")
# use random forest for imputation
meth[multicat] <- "cart"
# not imputed data
meth[c("subject_id", "survival_days", "event")] <- ""
# ───────────────────────────────────────────────────────────
# 3. predictorMatrix：
# ───────────────────────────────────────────────────────────
pred <- make.predictorMatrix(df_bf_impute)
pred[ , c("subject_id", "survival_days", "event",
          grep("^PO_", names(df_bf_impute), value = TRUE))] <- 0   # excluded from predictors
# ───────────────────────────────────────────────────────────
# 4. mice
# ───────────────────────────────────────────────────────────
set.seed(2025)
imp <- mice(df_bf_impute,
            m                  = 5,
            method             = meth,
            predictorMatrix    = pred,
            maxit              = 20,
            printFlag          = TRUE)







                    #### HANRUI's CODE END####









                    #### TRYING W MICE START--- IGNORE NO GOOD
### ___________________________________________________________________________

### impute survival days for subject_ids w/ survival_days = 90 and db_source == 
### "metavision"

### the imputed survival_days must be >= 90

### impute survival days for subject_ids w/ survival_days = 90 and db_source ==
### "metavision"
### the imputed survival_days must be >= 90
### Get completed dataset from previous MICE imputation
completed_data <- complete(imp, 1)


full_fixed_for_imputation = completed_data %>% 
  filter(survival_days >= 90) 


full_fixed_for_imputation %>% 
  filter(survival_days == 1460) %>% 
  count()


full_fixed_for_imputation %>% 
  filter(survival_days == 90) %>% 
  count()
### Identify subjects that need survival_days imputation
subjects_to_impute <- full_fixed_for_imputation %>%
  filter(survival_days == 90 & dbsource == "metavision") %>%
  pull(subject_id)
### Create dataset for survival imputation (set censored cases to NA)
surv_impute_data <- full_fixed_for_imputation
surv_impute_data$survival_days[surv_impute_data$subject_id %in% subjects_to_impute] <- NA
# Set up method for random forest imputation
meth1 <- make.method(surv_impute_data)
# Only use relevant variables as predictors for survival_days
surv_predictors <- c("age.at.admit", "gcs_total", "ethnicity",
                     "admission_type", "discharge_location", "event")
# Only impute survival_days using random forest
meth1["survival_days"] <- "pmm"  # Use random forest for survival_days
meth1[surv_predictors] <- "cart"
meth1[c("subject_id", "age.at.admit", "event", "gcs_total", "ethnicity",
        "admission_type", "discharge_location")] <- ""
# ───────────────────────────────────────────────────────────
# 3. predictorMatrix：
# ───────────────────────────────────────────────────────────
### Set up predictor matrix
pred1 <- make.predictorMatrix(surv_impute_data)
pred1[ , c("subject_id", "age.at.admit", "event", "gcs_total", "ethnicity",
           "admission_type", "discharge_location",
           grep("^PO_", names(surv_impute_data), value = TRUE))] <- 0   # excluded from predictors
# ───────────────────────────────────────────────────────────
# 4. mice
# ───────────────────────────────────────────────────────────
### Perform MICE with random forest
set.seed(2025)
surv_imp <- mice(surv_impute_data,
                 m               = 5,
                 method          = meth1,
                 predictorMatrix = pred1,
                 maxit           = 20,
                 printFlag       = TRUE)
surv_completed <- complete(surv_imp, 1)
### Ensure imputed values are >= 90
imputed_ids <- subjects_to_impute
for(id in imputed_ids) {
  current_val <- surv_completed$survival_days[surv_completed$subject_id == id]
  if(current_val < 90) {
    surv_completed$survival_days[surv_completed$subject_id == id] <- max(current_val, 90)
  }
}
### Final dataset with all imputations
df_final_imputed <- surv_completed
cat("Survival days imputation completed using random forest.\n")
if(length(subjects_to_impute) > 0) {
  imputed_vals <- df_final_imputed$survival_days[df_final_imputed$subject_id %in% subjects_to_impute]
  cat("Range of imputed survival_days:", min(imputed_vals), "to", max(imputed_vals), "\n")
}
View(df_final_imputed)
### Checking?
df_bf_impute %>%
  filter(survival_days == 90 & dbsource == "metavision") %>%
  select(subject_id, survival_days) %>% 
  View()
df_final_imputed %>%
  filter(dbsource == "metavision") %>%
  select(subject_id, survival_days) %>% 
  View()






              #### TRYING W MICE END--- IGNORE NO GOOD
### ___________________________________________________________________________























                           ### IT-MI demo:
                    ### BC MICE DOESNT WORK W THIS




### ___________________________________________________________________________

pacman::p_load(dplyr, tidyr, haven, 
               survival, ggplot2, ggfortify, 
               survminer, viridis, pscl, 
               extrafont, gridExtra, pscl,
               MASS, lubridate, lme4, lmtest,
               ggsurvfit, gdata, pseudo
)


dat <- complete(imp, 1)

dat_pseud <- dat %>%
  mutate(
    # Create sequential IDs
    ID = subject_id,  
    # Observed survival time
    Xi = survival_days,  
    # Create delta: 1 if event occurred, 0 if censored
    delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                   ifelse(survival_days == 1460,1,event))
  )


dat_pseud$ord_ID = 1:nrow(dat_pseud)


### patients censored @ 90 for metavision
tp.interest = 1460



my_dat = model.matrix(~comor.cancer+
                        drug.Parkinsons.Dementia+
                        language+
                        gcs_total+
                        age.at.admit+
                        proc.Cardiovasc+
                        gender+
                        religion, 
                      data = dat_pseud)







dat_pseud$po = pseudomean(dat_pseud$Xi, dat_pseud$delta, tmax=tp.interest)
baseline.mod = lm(dat_pseud$po~.-1, as.data.frame(my_dat))
# gets an unbiased estimate of min(Ti, tau) in the presence of censoring!!






summary(baseline.mod)
# now we will use this in the IT_MI algorithm that I pre-programmed for you.
source("C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/it_MI.R")






n = nrow(dat_pseud)
M = 5
ID = dat_pseud$ord_ID
Z = my_dat
  
  
  
X = dat_pseud$Xi
T_t = matrix(dat_pseud[, "Xi"], nrow = n, ncol = 1)
delta = matrix(dat_pseud[,"delta"], nrow = n, ncol = 1)
fitted.values = baseline.mod$coefficients





#  function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest)
# now to the it_MI file!!




MI5 = getMI(n, M, ID, Z, X, T_t, delta, fitted.values, tp.interest)
View(cbind(T_t, delta,MI5[[1]],MI5[[2]])) # view the first two imputed dataset))
# then you need to combine based on the SSIDs

# We have our imputed datasets in MI5, which is a list of matrices.
# Each matrix corresponds to an imputed dataset.

# combine each imputed dataset with the original data by making a column for 
# each imputed dataset's imputed survival days.

imputed_datasets <- lapply(MI5, function(mat) {
  data.frame(subject_id = dat_pseud$ID, imputed_survival_days = mat[,1])
})
# Combine all imputed datasets into one data frame
combined_imputed <- do.call(rbind, imputed_datasets)
# Add a column to indicate the imputation number
combined_imputed$imputation_number <- rep(1:M, each = nrow(dat_pseud))
# Now we can view the combined imputed datasets
View(combined_imputed)















