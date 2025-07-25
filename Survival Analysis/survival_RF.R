##### Ranger Survival Forest ####
pacman::p_load(survival, ranger, dplyr)
# taken from https://gist.github.com/thomasmooon/6eb87964ea663f4a7441cc2b2b730bd4


# mean imputation:
data = survival::cancer %>% mutate(status = status-1) %>% # 0 = censored, 1 = dead
  mutate(across(c(wt.loss, pat.karno, meal.cal, ph.ecog, ph.karno, age, sex, inst), ~ if_else(is.na(.), mean(., na.rm = T), .)),
         month.survival = time/30)
# some data contain missing values, for simplification I omit observations with NA's

test.ids = sample(nrow(data), size = 1/5*nrow(data))
train.ids = setdiff(1:nrow(data), test.ids)

test.data = data[test.ids, ]
train.data = data[train.ids , ]

rsf <- ranger(Surv(time = month.survival, event = status) ~ ., data = train.data)
pred.rsf = predict(rsf, data = test.data)
pred.rsf$unique.death.times
### gives you what you need to construct a kaplan-meier curve. 
# in this manner, we can get survival probability estimates at each time.
### Can get a survival curve.


# If we are interested in the E[T_i], we need to get an unbiased estimator
# of survival time first.
library(pseudo)
### must chose a survival time.
data = data %>% 
  mutate(pseudo.Ti = pseudomean(month.survival, status, 12))
mean(data$pseudo.Ti)
View(data)
test.data = data[test.ids, ]
train.data = data[train.ids , ]
pseudo.RF = ranger(pseudo.Ti ~ ., data = train.data %>% 
         dplyr::select(-time, -status))
View(pseudo.RF)
# outcomes are like that of a regression.


