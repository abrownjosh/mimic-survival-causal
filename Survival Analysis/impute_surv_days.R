### ROI before imputation 



                            #### HANRUI's CODE START####
rm(list = ls())


library(readr)
library(stringr)
library(mice)
library(dplyr)
library(glmnet)
library(Matrix)
library(survival)
data <- read.csv("C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/merged_v4.csv")
data <- data[,-1]

### MICE for multiple imputation after train-test splitting


# ───────────── 0 Preprocessing ───────────────────────────────
data <- data |>
  mutate(ethnicity = case_when(
    str_detect(ethnicity, regex("ASIAN"   , TRUE)) ~ "Asian",
    str_detect(ethnicity, regex("WHITE|MIDDLE", TRUE)) ~ "White",
    str_detect(ethnicity, regex("BLACK"   , TRUE)) ~ "Black or African American",
    str_detect(ethnicity, regex("OTHER|PATIENT|UNABLE|UNKNOWN", TRUE)) ~ "Unknown",
    str_detect(ethnicity, regex("AMERICAN", TRUE)) ~ "American Indian or Alaska Native",
    str_detect(ethnicity, regex("HISPANIC", TRUE)) ~ "Hispanic",
    str_detect(ethnicity, regex("MULTI"   , TRUE)) ~ "More than one race",
    TRUE ~ ethnicity)) |>
  mutate(across(where(is.character), factor))
# ───────────── 1 fix train / test ID  ────────────────
set.seed(2025)
all_id   <- unique(data$subject_id)
train_id <- sample(all_id, 0.60 * length(all_id))
test_id  <- setdiff(all_id, train_id)
train_raw <- data |> filter(subject_id %in% train_id)
test_raw  <- data |> filter(subject_id %in% test_id )
# ───────────── 2 mice settings ────────────────────────────
meth <- make.method(data)
meth[c("age.at.admit","gcs_total")] <- "pmm"
multicat <- c("admission_type","admission_location","discharge_location",
              "insurance","language","religion","marital_status","ethnicity")
meth[multicat] <- "cart"
meth[c("subject_id","survival_days","event")] <- ""
pred <- make.predictorMatrix(data)
pred[ , c("subject_id","survival_days","event",
          grep("^PO_", names(data), value = TRUE))] <- 0
# ───────────── 3 multiple imputation separately for train/test dataset ───────────────
imp_train <- mice(train_raw, m = 5, method = meth, predictorMatrix = pred,
                  maxit = 20, seed = 2025, printFlag = TRUE)






                    ##### QUINN DONT TOUCH TEST SET
imp_test  <- mice(test_raw , m = 5, method = meth, predictorMatrix = pred,
                  maxit = 20, seed = 2025, printFlag = TRUE)




train_list <- complete(imp_train, "all")   # list of 5 df
test_list  <- complete(imp_test , "all")

                    #### HANRUI's CODE END####





























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


dat <- complete(imp_train, 1)

#dat = train_list[[1]]


dat_pseud <- dat %>%
  mutate(
    # Create sequential IDs
    ID = subject_id,  
    # Observed survival time
    Xi = survival_days,  
    # Create delta: 1 if event occurred, 0 if censored
    delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                   ifelse(survival_days == 1460,1,event)),
    newReligion = ifelse(religion == "CATHOLIC", "CATHOLIC",
                         ifelse(religion == "NOT SPECIFIED", "NOT SPECIFIED",
                                ifelse(religion == "UNOBTAINABLE", "UNOBTAINABLE"
                                       , "OTHER"))),
    newLang = ifelse(language == "ENGL", "ENGL", "OTHER")
  )


dat_pseud$ord_ID = 1:nrow(dat_pseud)


### patients censored @ 90 for metavision
tp.interest = 1460



my_dat = model.matrix(~comor.cancer+
                        drug.Parkinsons.Dementia+
                        newLang+
                        gcs_total+
                        age.at.admit+
                        proc.Cardiovasc+
                        gender+
                        newReligion, 
                      data = dat_pseud)







dat_pseud$po = pseudomean(dat_pseud$Xi, dat_pseud$delta, tmax=tp.interest)
baseline.mod = lm(dat_pseud$po~.-1, as.data.frame(my_dat))
# gets an unbiased estimate of min(Ti, tau) in the presence of censoring!!






summary(baseline.mod)
# now we will use this in the IT_MI algorithm that I pre-programmed for you.
source("C:/Users/qqkoo/Downloads/it_MI.R")





n = nrow(dat_pseud)
M = 5
ID = dat_pseud$ord_ID
Z = my_dat
  
  
  
X = dat_pseud$Xi
T_t = matrix(dat_pseud[, "Xi"], nrow = n, ncol = 1)
delta = matrix(dat_pseud[,"delta"], nrow = n, ncol = 1)
fitted.values = baseline.mod$coefficients


#table(dat_pseud$newReligion)
#table(dat_pseud$newLang)

#  function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest)
# now to the it_MI file!!




MI5 = getMI(n, M, ID, Z, X, T_t, delta, fitted.values, tp.interest)
View(cbind(T_t, delta,MI5[[1]],MI5[[2]])) # view the first two imputed dataset))
# then you need to combine based on the SSIDs

# We have our imputed datasets in MI5, which is a list of matrices.
# Each matrix corresponds to an imputed dataset.qaq2






new1 = MI5[[1]]
new1 = data.frame(new1)

swag = cbind(dat,new1)

swag_fixed1 = swag %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                 ifelse(survival_days == 1460,1,event)))



swag_fixed1 = swag_fixed1 %>%
  mutate(surv_days_w_imps1 =ifelse(event==0,new1,survival_days)) %>% 
  dplyr::select(-new1) %>% 
  mutate(surv_days_w_imps1 = ifelse(surv_days_w_imps1 >1460, 1460, surv_days_w_imps1)) 





new2 = MI5[[2]]
new2 = data.frame(new2)

swag2 = cbind(swag_fixed1,new2)
swag_fixed2 = swag2 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                 ifelse(survival_days == 1460,1,event)))
swag_fixed2 = swag_fixed2 %>%
  mutate(surv_days_w_imps2 =ifelse(event==0,new2,survival_days)) %>% 
  dplyr::select(-new2) %>% 
  mutate(surv_days_w_imps2 = ifelse(surv_days_w_imps2 >1460, 1460, surv_days_w_imps2))





new3 = MI5[[3]]
new3 = data.frame(new3)
swag3 = cbind(swag_fixed2,new3)
swag_fixed3 = swag3 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                 ifelse(survival_days == 1460,1,event)))

swag_fixed3 = swag_fixed3 %>%
  mutate(surv_days_w_imps3 =ifelse(event==0,new3,survival_days)) %>% 
  dplyr::select(-new3) %>% 
  mutate(surv_days_w_imps3 = ifelse(surv_days_w_imps3 >1460, 1460, surv_days_w_imps3))




new4 = MI5[[4]]
new4 = data.frame(new4)
swag4 = cbind(swag_fixed3,new4)
swag_fixed4 = swag4 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                 ifelse(survival_days == 1460,1,event)))
swag_fixed4 = swag_fixed4 %>%
  mutate(surv_days_w_imps4 =ifelse(event==0,new4,survival_days)) %>% 
  dplyr::select(-new4) %>% 
  mutate(surv_days_w_imps4 = ifelse(surv_days_w_imps4 >1460, 1460, surv_days_w_imps4))

new5 = MI5[[5]]
new5 = data.frame(new5)
swag5 = cbind(swag_fixed4,new5)
swag_fixed5 = swag5 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                 ifelse(survival_days == 1460,1,event)))
swag_fixed5 = swag_fixed5 %>%
  mutate(surv_days_w_imps5 =ifelse(event==0,new5,survival_days)) %>% 
  dplyr::select(-new5) %>% 
  mutate(surv_days_w_imps5 = ifelse(surv_days_w_imps5 >1460, 1460, surv_days_w_imps5))
  
  
final_w_imps_1 = swag_fixed5 %>%
  group_by(subject_id) %>% 
  mutate(surv_days_w_mean_imps = sum(surv_days_w_imps1,
                                      surv_days_w_imps2,
                                      surv_days_w_imps3,
                                      surv_days_w_imps4,
                                      surv_days_w_imps5)/5) 
final_w_imps_1 = final_w_imps_1 %>%
  dplyr::select(-surv_days_w_imps1,
                -surv_days_w_imps2,
                -surv_days_w_imps3,
                -surv_days_w_imps4,
                -surv_days_w_imps5) %>%
  ungroup()

write.csv(final_w_imps_1, "C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/final_w_imps_1.csv", row.names = FALSE)
  


### ___________________________________________________________________________
# Now we will do the same thing for the second imputed dataset.
  
pacman::p_load(dplyr, tidyr, haven, 
               survival, ggplot2, ggfortify, 
               survminer, viridis, pscl, 
               extrafont, gridExtra, pscl,
               MASS, lubridate, lme4, lmtest,
               ggsurvfit, gdata, pseudo
)


dat <- complete(imp_train, 2)

#dat = train_list[[1]]


dat_pseud <- dat %>%
  mutate(
    # Create sequential IDs
    ID = subject_id,  
    # Observed survival time
    Xi = survival_days,  
    # Create delta: 1 if event occurred, 0 if censored
    delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                   ifelse(survival_days == 1460,1,event)),
    newReligion = ifelse(religion == "CATHOLIC", "CATHOLIC",
                         ifelse(religion == "NOT SPECIFIED", "NOT SPECIFIED",
                                ifelse(religion == "UNOBTAINABLE", "UNOBTAINABLE"
                                       , "OTHER"))),
    newLang = ifelse(language == "ENGL", "ENGL", "OTHER")
  )


dat_pseud$ord_ID = 1:nrow(dat_pseud)


### patients censored @ 90 for metavision
tp.interest = 1460



my_dat = model.matrix(~comor.cancer+
                        drug.Parkinsons.Dementia+
                        newLang+
                        gcs_total+
                        age.at.admit+
                        proc.Cardiovasc+
                        gender+
                        newReligion, 
                      data = dat_pseud)







dat_pseud$po = pseudomean(dat_pseud$Xi, dat_pseud$delta, tmax=tp.interest)
baseline.mod = lm(dat_pseud$po~.-1, as.data.frame(my_dat))
# gets an unbiased estimate of min(Ti, tau) in the presence of censoring!!






summary(baseline.mod)
# now we will use this in the IT_MI algorithm that I pre-programmed for you.
source("C:/Users/qqkoo/Downloads/it_MI.R")





n = nrow(dat_pseud)
M = 5
ID = dat_pseud$ord_ID
Z = my_dat



X = dat_pseud$Xi
T_t = matrix(dat_pseud[, "Xi"], nrow = n, ncol = 1)
delta = matrix(dat_pseud[,"delta"], nrow = n, ncol = 1)
fitted.values = baseline.mod$coefficients


#table(dat_pseud$newReligion)
#table(dat_pseud$newLang)

#  function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest)
# now to the it_MI file!!




MI5 = getMI(n, M, ID, Z, X, T_t, delta, fitted.values, tp.interest)
View(cbind(T_t, delta,MI5[[1]],MI5[[2]])) # view the first two imputed dataset))
# then you need to combine based on the SSIDs

# We have our imputed datasets in MI5, which is a list of matrices.
# Each matrix corresponds to an imputed dataset.qaq2






new1 = MI5[[1]]
new1 = data.frame(new1)

swag = cbind(dat,new1)

swag_fixed1 = swag %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))



swag_fixed1 = swag_fixed1 %>%
  mutate(surv_days_w_imps1 =ifelse(event==0,new1,survival_days)) %>% 
  dplyr::select(-new1) %>% 
  mutate(surv_days_w_imps1 = ifelse(surv_days_w_imps1 >1460, 1460, surv_days_w_imps1)) 





new2 = MI5[[2]]
new2 = data.frame(new2)

swag2 = cbind(swag_fixed1,new2)
swag_fixed2 = swag2 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed2 = swag_fixed2 %>%
  mutate(surv_days_w_imps2 =ifelse(event==0,new2,survival_days)) %>% 
  dplyr::select(-new2) %>% 
  mutate(surv_days_w_imps2 = ifelse(surv_days_w_imps2 >1460, 1460, surv_days_w_imps2))





new3 = MI5[[3]]
new3 = data.frame(new3)
swag3 = cbind(swag_fixed2,new3)
swag_fixed3 = swag3 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))

swag_fixed3 = swag_fixed3 %>%
  mutate(surv_days_w_imps3 =ifelse(event==0,new3,survival_days)) %>% 
  dplyr::select(-new3) %>% 
  mutate(surv_days_w_imps3 = ifelse(surv_days_w_imps3 >1460, 1460, surv_days_w_imps3))




new4 = MI5[[4]]
new4 = data.frame(new4)
swag4 = cbind(swag_fixed3,new4)
swag_fixed4 = swag4 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed4 = swag_fixed4 %>%
  mutate(surv_days_w_imps4 =ifelse(event==0,new4,survival_days)) %>% 
  dplyr::select(-new4) %>% 
  mutate(surv_days_w_imps4 = ifelse(surv_days_w_imps4 >1460, 1460, surv_days_w_imps4))

new5 = MI5[[5]]
new5 = data.frame(new5)
swag5 = cbind(swag_fixed4,new5)
swag_fixed5 = swag5 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed5 = swag_fixed5 %>%
  mutate(surv_days_w_imps5 =ifelse(event==0,new5,survival_days)) %>% 
  dplyr::select(-new5) %>% 
  mutate(surv_days_w_imps5 = ifelse(surv_days_w_imps5 >1460, 1460, surv_days_w_imps5))


final_w_imps_2 = swag_fixed5 %>%
  group_by(subject_id) %>% 
  mutate(surv_days_w_mean_imps = sum(surv_days_w_imps1,
                                     surv_days_w_imps2,
                                     surv_days_w_imps3,
                                     surv_days_w_imps4,
                                     surv_days_w_imps5)/5) 
final_w_imps_2 = final_w_imps_2 %>%
  dplyr::select(-surv_days_w_imps1,
                -surv_days_w_imps2,
                -surv_days_w_imps3,
                -surv_days_w_imps4,
                -surv_days_w_imps5) %>%
  ungroup()

write.csv(final_w_imps_2, "C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/final_w_imps_2.csv", row.names = FALSE)








### ___________________________________________________________________________
# Now we will do the same thing for the third imputed dataset.


pacman::p_load(dplyr, tidyr, haven, 
               survival, ggplot2, ggfortify, 
               survminer, viridis, pscl, 
               extrafont, gridExtra, pscl,
               MASS, lubridate, lme4, lmtest,
               ggsurvfit, gdata, pseudo
)


dat <- complete(imp_train, 3)

#dat = train_list[[1]]


dat_pseud <- dat %>%
  mutate(
    # Create sequential IDs
    ID = subject_id,  
    # Observed survival time
    Xi = survival_days,  
    # Create delta: 1 if event occurred, 0 if censored
    delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                   ifelse(survival_days == 1460,1,event)),
    newReligion = ifelse(religion == "CATHOLIC", "CATHOLIC",
                         ifelse(religion == "NOT SPECIFIED", "NOT SPECIFIED",
                                ifelse(religion == "UNOBTAINABLE", "UNOBTAINABLE"
                                       , "OTHER"))),
    newLang = ifelse(language == "ENGL", "ENGL", "OTHER")
  )


dat_pseud$ord_ID = 1:nrow(dat_pseud)


### patients censored @ 90 for metavision
tp.interest = 1460



my_dat = model.matrix(~comor.cancer+
                        drug.Parkinsons.Dementia+
                        newLang+
                        gcs_total+
                        age.at.admit+
                        proc.Cardiovasc+
                        gender+
                        newReligion, 
                      data = dat_pseud)







dat_pseud$po = pseudomean(dat_pseud$Xi, dat_pseud$delta, tmax=tp.interest)
baseline.mod = lm(dat_pseud$po~.-1, as.data.frame(my_dat))
# gets an unbiased estimate of min(Ti, tau) in the presence of censoring!!






summary(baseline.mod)
# now we will use this in the IT_MI algorithm that I pre-programmed for you.
source("C:/Users/qqkoo/Downloads/it_MI.R")





n = nrow(dat_pseud)
M = 5
ID = dat_pseud$ord_ID
Z = my_dat



X = dat_pseud$Xi
T_t = matrix(dat_pseud[, "Xi"], nrow = n, ncol = 1)
delta = matrix(dat_pseud[,"delta"], nrow = n, ncol = 1)
fitted.values = baseline.mod$coefficients


#table(dat_pseud$newReligion)
#table(dat_pseud$newLang)

#  function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest)
# now to the it_MI file!!




MI5 = getMI(n, M, ID, Z, X, T_t, delta, fitted.values, tp.interest)
View(cbind(T_t, delta,MI5[[1]],MI5[[2]])) # view the first two imputed dataset))
# then you need to combine based on the SSIDs

# We have our imputed datasets in MI5, which is a list of matrices.
# Each matrix corresponds to an imputed dataset.qaq2






new1 = MI5[[1]]
new1 = data.frame(new1)

swag = cbind(dat,new1)

swag_fixed1 = swag %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))



swag_fixed1 = swag_fixed1 %>%
  mutate(surv_days_w_imps1 =ifelse(event==0,new1,survival_days)) %>% 
  dplyr::select(-new1) %>% 
  mutate(surv_days_w_imps1 = ifelse(surv_days_w_imps1 >1460, 1460, surv_days_w_imps1)) 





new2 = MI5[[2]]
new2 = data.frame(new2)

swag2 = cbind(swag_fixed1,new2)
swag_fixed2 = swag2 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed2 = swag_fixed2 %>%
  mutate(surv_days_w_imps2 =ifelse(event==0,new2,survival_days)) %>% 
  dplyr::select(-new2) %>% 
  mutate(surv_days_w_imps2 = ifelse(surv_days_w_imps2 >1460, 1460, surv_days_w_imps2))





new3 = MI5[[3]]
new3 = data.frame(new3)
swag3 = cbind(swag_fixed2,new3)
swag_fixed3 = swag3 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))

swag_fixed3 = swag_fixed3 %>%
  mutate(surv_days_w_imps3 =ifelse(event==0,new3,survival_days)) %>% 
  dplyr::select(-new3) %>% 
  mutate(surv_days_w_imps3 = ifelse(surv_days_w_imps3 >1460, 1460, surv_days_w_imps3))




new4 = MI5[[4]]
new4 = data.frame(new4)
swag4 = cbind(swag_fixed3,new4)
swag_fixed4 = swag4 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed4 = swag_fixed4 %>%
  mutate(surv_days_w_imps4 =ifelse(event==0,new4,survival_days)) %>% 
  dplyr::select(-new4) %>% 
  mutate(surv_days_w_imps4 = ifelse(surv_days_w_imps4 >1460, 1460, surv_days_w_imps4))

new5 = MI5[[5]]
new5 = data.frame(new5)
swag5 = cbind(swag_fixed4,new5)
swag_fixed5 = swag5 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed5 = swag_fixed5 %>%
  mutate(surv_days_w_imps5 =ifelse(event==0,new5,survival_days)) %>% 
  dplyr::select(-new5) %>% 
  mutate(surv_days_w_imps5 = ifelse(surv_days_w_imps5 >1460, 1460, surv_days_w_imps5))


final_w_imps_3 = swag_fixed5 %>%
  group_by(subject_id) %>% 
  mutate(surv_days_w_mean_imps = sum(surv_days_w_imps1,
                                     surv_days_w_imps2,
                                     surv_days_w_imps3,
                                     surv_days_w_imps4,
                                     surv_days_w_imps5)/5) 
final_w_imps_3 = final_w_imps_3 %>%
  dplyr::select(-surv_days_w_imps1,
                -surv_days_w_imps2,
                -surv_days_w_imps3,
                -surv_days_w_imps4,
                -surv_days_w_imps5) %>%
  ungroup()

write.csv(final_w_imps_3, "C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/final_w_imps_3.csv", row.names = FALSE)

### ___________________________________________________________________________
# Now we will do the same thing for the fourth imputed dataset.

pacman::p_load(dplyr, tidyr, haven, 
               survival, ggplot2, ggfortify, 
               survminer, viridis, pscl, 
               extrafont, gridExtra, pscl,
               MASS, lubridate, lme4, lmtest,
               ggsurvfit, gdata, pseudo
)


dat <- complete(imp_train, 4)

#dat = train_list[[1]]


dat_pseud <- dat %>%
  mutate(
    # Create sequential IDs
    ID = subject_id,  
    # Observed survival time
    Xi = survival_days,  
    # Create delta: 1 if event occurred, 0 if censored
    delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                   ifelse(survival_days == 1460,1,event)),
    newReligion = ifelse(religion == "CATHOLIC", "CATHOLIC",
                         ifelse(religion == "NOT SPECIFIED", "NOT SPECIFIED",
                                ifelse(religion == "UNOBTAINABLE", "UNOBTAINABLE"
                                       , "OTHER"))),
    newLang = ifelse(language == "ENGL", "ENGL", "OTHER")
  )


dat_pseud$ord_ID = 1:nrow(dat_pseud)


### patients censored @ 90 for metavision
tp.interest = 1460



my_dat = model.matrix(~comor.cancer+
                        drug.Parkinsons.Dementia+
                        newLang+
                        gcs_total+
                        age.at.admit+
                        proc.Cardiovasc+
                        gender+
                        newReligion, 
                      data = dat_pseud)







dat_pseud$po = pseudomean(dat_pseud$Xi, dat_pseud$delta, tmax=tp.interest)
baseline.mod = lm(dat_pseud$po~.-1, as.data.frame(my_dat))
# gets an unbiased estimate of min(Ti, tau) in the presence of censoring!!






summary(baseline.mod)
# now we will use this in the IT_MI algorithm that I pre-programmed for you.
source("C:/Users/qqkoo/Downloads/it_MI.R")





n = nrow(dat_pseud)
M = 5
ID = dat_pseud$ord_ID
Z = my_dat



X = dat_pseud$Xi
T_t = matrix(dat_pseud[, "Xi"], nrow = n, ncol = 1)
delta = matrix(dat_pseud[,"delta"], nrow = n, ncol = 1)
fitted.values = baseline.mod$coefficients


#table(dat_pseud$newReligion)
#table(dat_pseud$newLang)

#  function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest)
# now to the it_MI file!!




MI5 = getMI(n, M, ID, Z, X, T_t, delta, fitted.values, tp.interest)
View(cbind(T_t, delta,MI5[[1]],MI5[[2]])) # view the first two imputed dataset))
# then you need to combine based on the SSIDs

# We have our imputed datasets in MI5, which is a list of matrices.
# Each matrix corresponds to an imputed dataset.qaq2






new1 = MI5[[1]]
new1 = data.frame(new1)

swag = cbind(dat,new1)

swag_fixed1 = swag %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))



swag_fixed1 = swag_fixed1 %>%
  mutate(surv_days_w_imps1 =ifelse(event==0,new1,survival_days)) %>% 
  dplyr::select(-new1) %>% 
  mutate(surv_days_w_imps1 = ifelse(surv_days_w_imps1 >1460, 1460, surv_days_w_imps1)) 





new2 = MI5[[2]]
new2 = data.frame(new2)

swag2 = cbind(swag_fixed1,new2)
swag_fixed2 = swag2 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed2 = swag_fixed2 %>%
  mutate(surv_days_w_imps2 =ifelse(event==0,new2,survival_days)) %>% 
  dplyr::select(-new2) %>% 
  mutate(surv_days_w_imps2 = ifelse(surv_days_w_imps2 >1460, 1460, surv_days_w_imps2))





new3 = MI5[[3]]
new3 = data.frame(new3)
swag3 = cbind(swag_fixed2,new3)
swag_fixed3 = swag3 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))

swag_fixed3 = swag_fixed3 %>%
  mutate(surv_days_w_imps3 =ifelse(event==0,new3,survival_days)) %>% 
  dplyr::select(-new3) %>% 
  mutate(surv_days_w_imps3 = ifelse(surv_days_w_imps3 >1460, 1460, surv_days_w_imps3))




new4 = MI5[[4]]
new4 = data.frame(new4)
swag4 = cbind(swag_fixed3,new4)
swag_fixed4 = swag4 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed4 = swag_fixed4 %>%
  mutate(surv_days_w_imps4 =ifelse(event==0,new4,survival_days)) %>% 
  dplyr::select(-new4) %>% 
  mutate(surv_days_w_imps4 = ifelse(surv_days_w_imps4 >1460, 1460, surv_days_w_imps4))

new5 = MI5[[5]]
new5 = data.frame(new5)
swag5 = cbind(swag_fixed4,new5)
swag_fixed5 = swag5 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed5 = swag_fixed5 %>%
  mutate(surv_days_w_imps5 =ifelse(event==0,new5,survival_days)) %>% 
  dplyr::select(-new5) %>% 
  mutate(surv_days_w_imps5 = ifelse(surv_days_w_imps5 >1460, 1460, surv_days_w_imps5))


final_w_imps_4 = swag_fixed5 %>%
  group_by(subject_id) %>% 
  mutate(surv_days_w_mean_imps = sum(surv_days_w_imps1,
                                     surv_days_w_imps2,
                                     surv_days_w_imps3,
                                     surv_days_w_imps4,
                                     surv_days_w_imps5)/5) 
final_w_imps_4 = final_w_imps_4 %>%
  dplyr::select(-surv_days_w_imps1,
                -surv_days_w_imps2,
                -surv_days_w_imps3,
                -surv_days_w_imps4,
                -surv_days_w_imps5) %>%
  ungroup()

write.csv(final_w_imps_4, "C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/final_w_imps_4.csv", row.names = FALSE)

### ___________________________________________________________________________
# Now we will do the same thing for the fifth imputed dataset.

pacman::p_load(dplyr, tidyr, haven, 
               survival, ggplot2, ggfortify, 
               survminer, viridis, pscl, 
               extrafont, gridExtra, pscl,
               MASS, lubridate, lme4, lmtest,
               ggsurvfit, gdata, pseudo
)


dat <- complete(imp_train, 5)

#dat = train_list[[1]]


dat_pseud <- dat %>%
  mutate(
    # Create sequential IDs
    ID = subject_id,  
    # Observed survival time
    Xi = survival_days,  
    # Create delta: 1 if event occurred, 0 if censored
    delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                   ifelse(survival_days == 1460,1,event)),
    newReligion = ifelse(religion == "CATHOLIC", "CATHOLIC",
                         ifelse(religion == "NOT SPECIFIED", "NOT SPECIFIED",
                                ifelse(religion == "UNOBTAINABLE", "UNOBTAINABLE"
                                       , "OTHER"))),
    newLang = ifelse(language == "ENGL", "ENGL", "OTHER")
  )


dat_pseud$ord_ID = 1:nrow(dat_pseud)


### patients censored @ 90 for metavision
tp.interest = 1460



my_dat = model.matrix(~comor.cancer+
                        drug.Parkinsons.Dementia+
                        newLang+
                        gcs_total+
                        age.at.admit+
                        proc.Cardiovasc+
                        gender+
                        newReligion, 
                      data = dat_pseud)







dat_pseud$po = pseudomean(dat_pseud$Xi, dat_pseud$delta, tmax=tp.interest)
baseline.mod = lm(dat_pseud$po~.-1, as.data.frame(my_dat))
# gets an unbiased estimate of min(Ti, tau) in the presence of censoring!!






summary(baseline.mod)
# now we will use this in the IT_MI algorithm that I pre-programmed for you.
source("C:/Users/qqkoo/Downloads/it_MI.R")





n = nrow(dat_pseud)
M = 5
ID = dat_pseud$ord_ID
Z = my_dat



X = dat_pseud$Xi
T_t = matrix(dat_pseud[, "Xi"], nrow = n, ncol = 1)
delta = matrix(dat_pseud[,"delta"], nrow = n, ncol = 1)
fitted.values = baseline.mod$coefficients


#table(dat_pseud$newReligion)
#table(dat_pseud$newLang)

#  function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest)
# now to the it_MI file!!




MI5 = getMI(n, M, ID, Z, X, T_t, delta, fitted.values, tp.interest)
View(cbind(T_t, delta,MI5[[1]],MI5[[2]])) # view the first two imputed dataset))
# then you need to combine based on the SSIDs

# We have our imputed datasets in MI5, which is a list of matrices.
# Each matrix corresponds to an imputed dataset.qaq2






new1 = MI5[[1]]
new1 = data.frame(new1)

swag = cbind(dat,new1)

swag_fixed1 = swag %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))



swag_fixed1 = swag_fixed1 %>%
  mutate(surv_days_w_imps1 =ifelse(event==0,new1,survival_days)) %>% 
  dplyr::select(-new1) %>% 
  mutate(surv_days_w_imps1 = ifelse(surv_days_w_imps1 >1460, 1460, surv_days_w_imps1)) 





new2 = MI5[[2]]
new2 = data.frame(new2)

swag2 = cbind(swag_fixed1,new2)
swag_fixed2 = swag2 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed2 = swag_fixed2 %>%
  mutate(surv_days_w_imps2 =ifelse(event==0,new2,survival_days)) %>% 
  dplyr::select(-new2) %>% 
  mutate(surv_days_w_imps2 = ifelse(surv_days_w_imps2 >1460, 1460, surv_days_w_imps2))





new3 = MI5[[3]]
new3 = data.frame(new3)
swag3 = cbind(swag_fixed2,new3)
swag_fixed3 = swag3 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))

swag_fixed3 = swag_fixed3 %>%
  mutate(surv_days_w_imps3 =ifelse(event==0,new3,survival_days)) %>% 
  dplyr::select(-new3) %>% 
  mutate(surv_days_w_imps3 = ifelse(surv_days_w_imps3 >1460, 1460, surv_days_w_imps3))




new4 = MI5[[4]]
new4 = data.frame(new4)
swag4 = cbind(swag_fixed3,new4)
swag_fixed4 = swag4 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed4 = swag_fixed4 %>%
  mutate(surv_days_w_imps4 =ifelse(event==0,new4,survival_days)) %>% 
  dplyr::select(-new4) %>% 
  mutate(surv_days_w_imps4 = ifelse(surv_days_w_imps4 >1460, 1460, surv_days_w_imps4))

new5 = MI5[[5]]
new5 = data.frame(new5)
swag5 = cbind(swag_fixed4,new5)
swag_fixed5 = swag5 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed5 = swag_fixed5 %>%
  mutate(surv_days_w_imps5 =ifelse(event==0,new5,survival_days)) %>% 
  dplyr::select(-new5) %>% 
  mutate(surv_days_w_imps5 = ifelse(surv_days_w_imps5 >1460, 1460, surv_days_w_imps5))


final_w_imps_5 = swag_fixed5 %>%
  group_by(subject_id) %>% 
  mutate(surv_days_w_mean_imps = sum(surv_days_w_imps1,
                                     surv_days_w_imps2,
                                     surv_days_w_imps3,
                                     surv_days_w_imps4,
                                     surv_days_w_imps5)/5) 
final_w_imps_5 = final_w_imps_5 %>%
  dplyr::select(-surv_days_w_imps1,
                -surv_days_w_imps2,
                -surv_days_w_imps3,
                -surv_days_w_imps4,
                -surv_days_w_imps5) %>%
  ungroup()

write.csv(final_w_imps_5, "C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/final_w_imps_5.csv", row.names = FALSE)











#### NOW DO FOR THE TEST SET




### ___________________________________________________________________________
pacman::p_load(dplyr, tidyr, haven, 
               survival, ggplot2, ggfortify, 
               survminer, viridis, pscl, 
               extrafont, gridExtra, pscl,
               MASS, lubridate, lme4, lmtest,
               ggsurvfit, gdata, pseudo
)


dat <- complete(imp_test, 1)

#dat = train_list[[1]]


dat_pseud <- dat %>%
  mutate(
    # Create sequential IDs
    ID = subject_id,  
    # Observed survival time
    Xi = survival_days,  
    # Create delta: 1 if event occurred, 0 if censored
    delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                   ifelse(survival_days == 1460,1,event)),
    newReligion = ifelse(religion == "CATHOLIC", "CATHOLIC",
                         ifelse(religion == "NOT SPECIFIED", "NOT SPECIFIED",
                                ifelse(religion == "UNOBTAINABLE", "UNOBTAINABLE"
                                       , "OTHER"))),
    newLang = ifelse(language == "ENGL", "ENGL", "OTHER")
  )


dat_pseud$ord_ID = 1:nrow(dat_pseud)


### patients censored @ 90 for metavision
tp.interest = 1460



my_dat = model.matrix(~comor.cancer+
                        drug.Parkinsons.Dementia+
                        newLang+
                        gcs_total+
                        age.at.admit+
                        proc.Cardiovasc+
                        gender+
                        newReligion, 
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


#table(dat_pseud$newReligion)
#table(dat_pseud$newLang)

#  function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest)
# now to the it_MI file!!




MI5 = getMI(n, M, ID, Z, X, T_t, delta, fitted.values, tp.interest)
View(cbind(T_t, delta,MI5[[1]],MI5[[2]])) # view the first two imputed dataset))
# then you need to combine based on the SSIDs

# We have our imputed datasets in MI5, which is a list of matrices.
# Each matrix corresponds to an imputed dataset.qaq2






new1 = MI5[[1]]
new1 = data.frame(new1)

swag = cbind(dat,new1)

swag_fixed1 = swag %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))



swag_fixed1 = swag_fixed1 %>%
  mutate(surv_days_w_imps1 =ifelse(event==0,new1,survival_days)) %>% 
  dplyr::select(-new1) %>% 
  mutate(surv_days_w_imps1 = ifelse(surv_days_w_imps1 >1460, 1460, surv_days_w_imps1)) 





new2 = MI5[[2]]
new2 = data.frame(new2)

swag2 = cbind(swag_fixed1,new2)
swag_fixed2 = swag2 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed2 = swag_fixed2 %>%
  mutate(surv_days_w_imps2 =ifelse(event==0,new2,survival_days)) %>% 
  dplyr::select(-new2) %>% 
  mutate(surv_days_w_imps2 = ifelse(surv_days_w_imps2 >1460, 1460, surv_days_w_imps2))





new3 = MI5[[3]]
new3 = data.frame(new3)
swag3 = cbind(swag_fixed2,new3)
swag_fixed3 = swag3 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))

swag_fixed3 = swag_fixed3 %>%
  mutate(surv_days_w_imps3 =ifelse(event==0,new3,survival_days)) %>% 
  dplyr::select(-new3) %>% 
  mutate(surv_days_w_imps3 = ifelse(surv_days_w_imps3 >1460, 1460, surv_days_w_imps3))




new4 = MI5[[4]]
new4 = data.frame(new4)
swag4 = cbind(swag_fixed3,new4)
swag_fixed4 = swag4 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed4 = swag_fixed4 %>%
  mutate(surv_days_w_imps4 =ifelse(event==0,new4,survival_days)) %>% 
  dplyr::select(-new4) %>% 
  mutate(surv_days_w_imps4 = ifelse(surv_days_w_imps4 >1460, 1460, surv_days_w_imps4))

new5 = MI5[[5]]
new5 = data.frame(new5)
swag5 = cbind(swag_fixed4,new5)
swag_fixed5 = swag5 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed5 = swag_fixed5 %>%
  mutate(surv_days_w_imps5 =ifelse(event==0,new5,survival_days)) %>% 
  dplyr::select(-new5) %>% 
  mutate(surv_days_w_imps5 = ifelse(surv_days_w_imps5 >1460, 1460, surv_days_w_imps5))


final_w_imps_1 = swag_fixed5 %>%
  group_by(subject_id) %>% 
  mutate(surv_days_w_mean_imps = sum(surv_days_w_imps1,
                                     surv_days_w_imps2,
                                     surv_days_w_imps3,
                                     surv_days_w_imps4,
                                     surv_days_w_imps5)/5) 
final_w_imps_1_test = final_w_imps_1 %>%
  dplyr::select(-surv_days_w_imps1,
                -surv_days_w_imps2,
                -surv_days_w_imps3,
                -surv_days_w_imps4,
                -surv_days_w_imps5) %>%
  ungroup()

write.csv(final_w_imps_1_test, "C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/final_w_imps_1_test.csv", row.names = FALSE)



### ___________________________________________________________________________
# Now we will do the same thing for the second imputed dataset.

pacman::p_load(dplyr, tidyr, haven, 
               survival, ggplot2, ggfortify, 
               survminer, viridis, pscl, 
               extrafont, gridExtra, pscl,
               MASS, lubridate, lme4, lmtest,
               ggsurvfit, gdata, pseudo
)


dat <- complete(imp_test, 2)

#dat = train_list[[1]]


dat_pseud <- dat %>%
  mutate(
    # Create sequential IDs
    ID = subject_id,  
    # Observed survival time
    Xi = survival_days,  
    # Create delta: 1 if event occurred, 0 if censored
    delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                   ifelse(survival_days == 1460,1,event)),
    newReligion = ifelse(religion == "CATHOLIC", "CATHOLIC",
                         ifelse(religion == "NOT SPECIFIED", "NOT SPECIFIED",
                                ifelse(religion == "UNOBTAINABLE", "UNOBTAINABLE"
                                       , "OTHER"))),
    newLang = ifelse(language == "ENGL", "ENGL", "OTHER")
  )


dat_pseud$ord_ID = 1:nrow(dat_pseud)


### patients censored @ 90 for metavision
tp.interest = 1460



my_dat = model.matrix(~comor.cancer+
                        drug.Parkinsons.Dementia+
                        newLang+
                        gcs_total+
                        age.at.admit+
                        proc.Cardiovasc+
                        gender+
                        newReligion, 
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


#table(dat_pseud$newReligion)
#table(dat_pseud$newLang)

#  function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest)
# now to the it_MI file!!




MI5 = getMI(n, M, ID, Z, X, T_t, delta, fitted.values, tp.interest)
View(cbind(T_t, delta,MI5[[1]],MI5[[2]])) # view the first two imputed dataset))
# then you need to combine based on the SSIDs

# We have our imputed datasets in MI5, which is a list of matrices.
# Each matrix corresponds to an imputed dataset.qaq2






new1 = MI5[[1]]
new1 = data.frame(new1)

swag = cbind(dat,new1)

swag_fixed1 = swag %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))



swag_fixed1 = swag_fixed1 %>%
  mutate(surv_days_w_imps1 =ifelse(event==0,new1,survival_days)) %>% 
  dplyr::select(-new1) %>% 
  mutate(surv_days_w_imps1 = ifelse(surv_days_w_imps1 >1460, 1460, surv_days_w_imps1)) 





new2 = MI5[[2]]
new2 = data.frame(new2)

swag2 = cbind(swag_fixed1,new2)
swag_fixed2 = swag2 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed2 = swag_fixed2 %>%
  mutate(surv_days_w_imps2 =ifelse(event==0,new2,survival_days)) %>% 
  dplyr::select(-new2) %>% 
  mutate(surv_days_w_imps2 = ifelse(surv_days_w_imps2 >1460, 1460, surv_days_w_imps2))





new3 = MI5[[3]]
new3 = data.frame(new3)
swag3 = cbind(swag_fixed2,new3)
swag_fixed3 = swag3 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))

swag_fixed3 = swag_fixed3 %>%
  mutate(surv_days_w_imps3 =ifelse(event==0,new3,survival_days)) %>% 
  dplyr::select(-new3) %>% 
  mutate(surv_days_w_imps3 = ifelse(surv_days_w_imps3 >1460, 1460, surv_days_w_imps3))




new4 = MI5[[4]]
new4 = data.frame(new4)
swag4 = cbind(swag_fixed3,new4)
swag_fixed4 = swag4 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed4 = swag_fixed4 %>%
  mutate(surv_days_w_imps4 =ifelse(event==0,new4,survival_days)) %>% 
  dplyr::select(-new4) %>% 
  mutate(surv_days_w_imps4 = ifelse(surv_days_w_imps4 >1460, 1460, surv_days_w_imps4))

new5 = MI5[[5]]
new5 = data.frame(new5)
swag5 = cbind(swag_fixed4,new5)
swag_fixed5 = swag5 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed5 = swag_fixed5 %>%
  mutate(surv_days_w_imps5 =ifelse(event==0,new5,survival_days)) %>% 
  dplyr::select(-new5) %>% 
  mutate(surv_days_w_imps5 = ifelse(surv_days_w_imps5 >1460, 1460, surv_days_w_imps5))


final_w_imps_2 = swag_fixed5 %>%
  group_by(subject_id) %>% 
  mutate(surv_days_w_mean_imps = sum(surv_days_w_imps1,
                                     surv_days_w_imps2,
                                     surv_days_w_imps3,
                                     surv_days_w_imps4,
                                     surv_days_w_imps5)/5) 
final_w_imps_2_test = final_w_imps_2 %>%
  dplyr::select(-surv_days_w_imps1,
                -surv_days_w_imps2,
                -surv_days_w_imps3,
                -surv_days_w_imps4,
                -surv_days_w_imps5) %>%
  ungroup()

write.csv(final_w_imps_2_test, "C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/final_w_imps_2_test.csv", row.names = FALSE)








### ___________________________________________________________________________
# Now we will do the same thing for the third imputed dataset.


pacman::p_load(dplyr, tidyr, haven, 
               survival, ggplot2, ggfortify, 
               survminer, viridis, pscl, 
               extrafont, gridExtra, pscl,
               MASS, lubridate, lme4, lmtest,
               ggsurvfit, gdata, pseudo
)


dat <- complete(imp_test, 3)

#dat = train_list[[1]]


dat_pseud <- dat %>%
  mutate(
    # Create sequential IDs
    ID = subject_id,  
    # Observed survival time
    Xi = survival_days,  
    # Create delta: 1 if event occurred, 0 if censored
    delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                   ifelse(survival_days == 1460,1,event)),
    newReligion = ifelse(religion == "CATHOLIC", "CATHOLIC",
                         ifelse(religion == "NOT SPECIFIED", "NOT SPECIFIED",
                                ifelse(religion == "UNOBTAINABLE", "UNOBTAINABLE"
                                       , "OTHER"))),
    newLang = ifelse(language == "ENGL", "ENGL", "OTHER")
  )


dat_pseud$ord_ID = 1:nrow(dat_pseud)


### patients censored @ 90 for metavision
tp.interest = 1460



my_dat = model.matrix(~comor.cancer+
                        drug.Parkinsons.Dementia+
                        newLang+
                        gcs_total+
                        age.at.admit+
                        proc.Cardiovasc+
                        gender+
                        newReligion, 
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


#table(dat_pseud$newReligion)
#table(dat_pseud$newLang)

#  function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest)
# now to the it_MI file!!




MI5 = getMI(n, M, ID, Z, X, T_t, delta, fitted.values, tp.interest)
View(cbind(T_t, delta,MI5[[1]],MI5[[2]])) # view the first two imputed dataset))
# then you need to combine based on the SSIDs

# We have our imputed datasets in MI5, which is a list of matrices.
# Each matrix corresponds to an imputed dataset.qaq2






new1 = MI5[[1]]
new1 = data.frame(new1)

swag = cbind(dat,new1)

swag_fixed1 = swag %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))



swag_fixed1 = swag_fixed1 %>%
  mutate(surv_days_w_imps1 =ifelse(event==0,new1,survival_days)) %>% 
  dplyr::select(-new1) %>% 
  mutate(surv_days_w_imps1 = ifelse(surv_days_w_imps1 >1460, 1460, surv_days_w_imps1)) 





new2 = MI5[[2]]
new2 = data.frame(new2)

swag2 = cbind(swag_fixed1,new2)
swag_fixed2 = swag2 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed2 = swag_fixed2 %>%
  mutate(surv_days_w_imps2 =ifelse(event==0,new2,survival_days)) %>% 
  dplyr::select(-new2) %>% 
  mutate(surv_days_w_imps2 = ifelse(surv_days_w_imps2 >1460, 1460, surv_days_w_imps2))





new3 = MI5[[3]]
new3 = data.frame(new3)
swag3 = cbind(swag_fixed2,new3)
swag_fixed3 = swag3 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))

swag_fixed3 = swag_fixed3 %>%
  mutate(surv_days_w_imps3 =ifelse(event==0,new3,survival_days)) %>% 
  dplyr::select(-new3) %>% 
  mutate(surv_days_w_imps3 = ifelse(surv_days_w_imps3 >1460, 1460, surv_days_w_imps3))




new4 = MI5[[4]]
new4 = data.frame(new4)
swag4 = cbind(swag_fixed3,new4)
swag_fixed4 = swag4 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed4 = swag_fixed4 %>%
  mutate(surv_days_w_imps4 =ifelse(event==0,new4,survival_days)) %>% 
  dplyr::select(-new4) %>% 
  mutate(surv_days_w_imps4 = ifelse(surv_days_w_imps4 >1460, 1460, surv_days_w_imps4))

new5 = MI5[[5]]
new5 = data.frame(new5)
swag5 = cbind(swag_fixed4,new5)
swag_fixed5 = swag5 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed5 = swag_fixed5 %>%
  mutate(surv_days_w_imps5 =ifelse(event==0,new5,survival_days)) %>% 
  dplyr::select(-new5) %>% 
  mutate(surv_days_w_imps5 = ifelse(surv_days_w_imps5 >1460, 1460, surv_days_w_imps5))


final_w_imps_3 = swag_fixed5 %>%
  group_by(subject_id) %>% 
  mutate(surv_days_w_mean_imps = sum(surv_days_w_imps1,
                                     surv_days_w_imps2,
                                     surv_days_w_imps3,
                                     surv_days_w_imps4,
                                     surv_days_w_imps5)/5) 
final_w_imps_3_test = final_w_imps_3 %>%
  dplyr::select(-surv_days_w_imps1,
                -surv_days_w_imps2,
                -surv_days_w_imps3,
                -surv_days_w_imps4,
                -surv_days_w_imps5) %>%
  ungroup()

write.csv(final_w_imps_3_test, "C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/final_w_imps_3_test.csv", row.names = FALSE)

### ___________________________________________________________________________
# Now we will do the same thing for the fourth imputed dataset.

pacman::p_load(dplyr, tidyr, haven, 
               survival, ggplot2, ggfortify, 
               survminer, viridis, pscl, 
               extrafont, gridExtra, pscl,
               MASS, lubridate, lme4, lmtest,
               ggsurvfit, gdata, pseudo
)


dat <- complete(imp_test, 4)

#dat = train_list[[1]]


dat_pseud <- dat %>%
  mutate(
    # Create sequential IDs
    ID = subject_id,  
    # Observed survival time
    Xi = survival_days,  
    # Create delta: 1 if event occurred, 0 if censored
    delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                   ifelse(survival_days == 1460,1,event)),
    newReligion = ifelse(religion == "CATHOLIC", "CATHOLIC",
                         ifelse(religion == "NOT SPECIFIED", "NOT SPECIFIED",
                                ifelse(religion == "UNOBTAINABLE", "UNOBTAINABLE"
                                       , "OTHER"))),
    newLang = ifelse(language == "ENGL", "ENGL", "OTHER")
  )


dat_pseud$ord_ID = 1:nrow(dat_pseud)


### patients censored @ 90 for metavision
tp.interest = 1460



my_dat = model.matrix(~comor.cancer+
                        drug.Parkinsons.Dementia+
                        newLang+
                        gcs_total+
                        age.at.admit+
                        proc.Cardiovasc+
                        gender+
                        newReligion, 
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


#table(dat_pseud$newReligion)
#table(dat_pseud$newLang)

#  function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest)
# now to the it_MI file!!




MI5 = getMI(n, M, ID, Z, X, T_t, delta, fitted.values, tp.interest)
View(cbind(T_t, delta,MI5[[1]],MI5[[2]])) # view the first two imputed dataset))
# then you need to combine based on the SSIDs

# We have our imputed datasets in MI5, which is a list of matrices.
# Each matrix corresponds to an imputed dataset.qaq2






new1 = MI5[[1]]
new1 = data.frame(new1)

swag = cbind(dat,new1)

swag_fixed1 = swag %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))



swag_fixed1 = swag_fixed1 %>%
  mutate(surv_days_w_imps1 =ifelse(event==0,new1,survival_days)) %>% 
  dplyr::select(-new1) %>% 
  mutate(surv_days_w_imps1 = ifelse(surv_days_w_imps1 >1460, 1460, surv_days_w_imps1)) 





new2 = MI5[[2]]
new2 = data.frame(new2)

swag2 = cbind(swag_fixed1,new2)
swag_fixed2 = swag2 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed2 = swag_fixed2 %>%
  mutate(surv_days_w_imps2 =ifelse(event==0,new2,survival_days)) %>% 
  dplyr::select(-new2) %>% 
  mutate(surv_days_w_imps2 = ifelse(surv_days_w_imps2 >1460, 1460, surv_days_w_imps2))





new3 = MI5[[3]]
new3 = data.frame(new3)
swag3 = cbind(swag_fixed2,new3)
swag_fixed3 = swag3 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))

swag_fixed3 = swag_fixed3 %>%
  mutate(surv_days_w_imps3 =ifelse(event==0,new3,survival_days)) %>% 
  dplyr::select(-new3) %>% 
  mutate(surv_days_w_imps3 = ifelse(surv_days_w_imps3 >1460, 1460, surv_days_w_imps3))




new4 = MI5[[4]]
new4 = data.frame(new4)
swag4 = cbind(swag_fixed3,new4)
swag_fixed4 = swag4 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed4 = swag_fixed4 %>%
  mutate(surv_days_w_imps4 =ifelse(event==0,new4,survival_days)) %>% 
  dplyr::select(-new4) %>% 
  mutate(surv_days_w_imps4 = ifelse(surv_days_w_imps4 >1460, 1460, surv_days_w_imps4))

new5 = MI5[[5]]
new5 = data.frame(new5)
swag5 = cbind(swag_fixed4,new5)
swag_fixed5 = swag5 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed5 = swag_fixed5 %>%
  mutate(surv_days_w_imps5 =ifelse(event==0,new5,survival_days)) %>% 
  dplyr::select(-new5) %>% 
  mutate(surv_days_w_imps5 = ifelse(surv_days_w_imps5 >1460, 1460, surv_days_w_imps5))


final_w_imps_4 = swag_fixed5 %>%
  group_by(subject_id) %>% 
  mutate(surv_days_w_mean_imps = sum(surv_days_w_imps1,
                                     surv_days_w_imps2,
                                     surv_days_w_imps3,
                                     surv_days_w_imps4,
                                     surv_days_w_imps5)/5) 
final_w_imps_4_test = final_w_imps_4 %>%
  dplyr::select(-surv_days_w_imps1,
                -surv_days_w_imps2,
                -surv_days_w_imps3,
                -surv_days_w_imps4,
                -surv_days_w_imps5) %>%
  ungroup()

write.csv(final_w_imps_4_test, "C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/final_w_imps_4_test.csv", row.names = FALSE)

### ___________________________________________________________________________
# Now we will do the same thing for the fifth imputed dataset.

pacman::p_load(dplyr, tidyr, haven, 
               survival, ggplot2, ggfortify, 
               survminer, viridis, pscl, 
               extrafont, gridExtra, pscl,
               MASS, lubridate, lme4, lmtest,
               ggsurvfit, gdata, pseudo
)


dat <- complete(imp_test, 5)

#dat = train_list[[1]]


dat_pseud <- dat %>%
  mutate(
    # Create sequential IDs
    ID = subject_id,  
    # Observed survival time
    Xi = survival_days,  
    # Create delta: 1 if event occurred, 0 if censored
    delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                   ifelse(survival_days == 1460,1,event)),
    newReligion = ifelse(religion == "CATHOLIC", "CATHOLIC",
                         ifelse(religion == "NOT SPECIFIED", "NOT SPECIFIED",
                                ifelse(religion == "UNOBTAINABLE", "UNOBTAINABLE"
                                       , "OTHER"))),
    newLang = ifelse(language == "ENGL", "ENGL", "OTHER")
  )


dat_pseud$ord_ID = 1:nrow(dat_pseud)


### patients censored @ 90 for metavision
tp.interest = 1460



my_dat = model.matrix(~comor.cancer+
                        drug.Parkinsons.Dementia+
                        newLang+
                        gcs_total+
                        age.at.admit+
                        proc.Cardiovasc+
                        gender+
                        newReligion, 
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


#table(dat_pseud$newReligion)
#table(dat_pseud$newLang)

#  function(n,M,ID,Z,X,T_t,delta, fitted.values, tp.interest)
# now to the it_MI file!!




MI5 = getMI(n, M, ID, Z, X, T_t, delta, fitted.values, tp.interest)
View(cbind(T_t, delta,MI5[[1]],MI5[[2]])) # view the first two imputed dataset))
# then you need to combine based on the SSIDs

# We have our imputed datasets in MI5, which is a list of matrices.
# Each matrix corresponds to an imputed dataset.qaq2






new1 = MI5[[1]]
new1 = data.frame(new1)

swag = cbind(dat,new1)

swag_fixed1 = swag %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))



swag_fixed1 = swag_fixed1 %>%
  mutate(surv_days_w_imps1 =ifelse(event==0,new1,survival_days)) %>% 
  dplyr::select(-new1) %>% 
  mutate(surv_days_w_imps1 = ifelse(surv_days_w_imps1 >1460, 1460, surv_days_w_imps1)) 





new2 = MI5[[2]]
new2 = data.frame(new2)

swag2 = cbind(swag_fixed1,new2)
swag_fixed2 = swag2 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed2 = swag_fixed2 %>%
  mutate(surv_days_w_imps2 =ifelse(event==0,new2,survival_days)) %>% 
  dplyr::select(-new2) %>% 
  mutate(surv_days_w_imps2 = ifelse(surv_days_w_imps2 >1460, 1460, surv_days_w_imps2))





new3 = MI5[[3]]
new3 = data.frame(new3)
swag3 = cbind(swag_fixed2,new3)
swag_fixed3 = swag3 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))

swag_fixed3 = swag_fixed3 %>%
  mutate(surv_days_w_imps3 =ifelse(event==0,new3,survival_days)) %>% 
  dplyr::select(-new3) %>% 
  mutate(surv_days_w_imps3 = ifelse(surv_days_w_imps3 >1460, 1460, surv_days_w_imps3))




new4 = MI5[[4]]
new4 = data.frame(new4)
swag4 = cbind(swag_fixed3,new4)
swag_fixed4 = swag4 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed4 = swag_fixed4 %>%
  mutate(surv_days_w_imps4 =ifelse(event==0,new4,survival_days)) %>% 
  dplyr::select(-new4) %>% 
  mutate(surv_days_w_imps4 = ifelse(surv_days_w_imps4 >1460, 1460, surv_days_w_imps4))

new5 = MI5[[5]]
new5 = data.frame(new5)
swag5 = cbind(swag_fixed4,new5)
swag_fixed5 = swag5 %>%
  mutate(delta = ifelse(event == 0 & dbsource == "metavision", 0, 
                        ifelse(survival_days == 1460,1,event)))
swag_fixed5 = swag_fixed5 %>%
  mutate(surv_days_w_imps5 =ifelse(event==0,new5,survival_days)) %>% 
  dplyr::select(-new5) %>% 
  mutate(surv_days_w_imps5 = ifelse(surv_days_w_imps5 >1460, 1460, surv_days_w_imps5))


final_w_imps_5 = swag_fixed5 %>%
  group_by(subject_id) %>% 
  mutate(surv_days_w_mean_imps = sum(surv_days_w_imps1,
                                     surv_days_w_imps2,
                                     surv_days_w_imps3,
                                     surv_days_w_imps4,
                                     surv_days_w_imps5)/5) 
final_w_imps_5_test = final_w_imps_5 %>%
  dplyr::select(-surv_days_w_imps1,
                -surv_days_w_imps2,
                -surv_days_w_imps3,
                -surv_days_w_imps4,
                -surv_days_w_imps5) %>%
  ungroup()

write.csv(final_w_imps_5_test, "C:/Users/qqkoo/OneDrive/Desktop/BDSI/ML_SOI_FINAL/final_w_imps_5_test.csv", row.names = FALSE)








