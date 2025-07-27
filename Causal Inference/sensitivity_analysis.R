knitr::opts_chunk$set(echo = TRUE)

# Clear the environment and load packages
rm(list = ls())
pacman::p_load(tidyverse, sensemakr)

# Set directory
mydir = "~/BDSI/research_group/data"
setwd(mydir)


# Data Definition, Reading, and Cleaning 
# _____________________________________________________________________________

chart.data = read.csv("imputed6.csv")

chart.data = chart.data %>% 
  rename(age = age.at.admit, gcs = gcs.initial)

chart.min = chart.data %>%
  filter(intervention.group != "Craniotomy/Craniectomy")

# Rename values to 0 and 1
chart.min = chart.min %>%
  mutate(intervention.group = if_else(intervention.group == 
                                        "Minimally Invasive Surgery", 1, 0))

chart.cranio = chart.data %>% 
  filter(intervention.group != "Minimally Invasive Surgery")

# Rename values to 0 and 1
chart.cranio = chart.cranio %>%
  mutate(intervention.group = if_else(intervention.group == 
                                        "Craniotomy/Craniectomy", 1, 0))


# Minimally Invasive Surgery Sensitivity Analysis
# _____________________________________________________________________________

# Linear regression model
lin_min = lm(died.before.90.days ~ intervention.group + age + 
                 gcs + wbc + platelet + hematocrit + glucose + 
                 creatinine + sodium + potassium + diabetes + anemia + 
                 liver.disease + depressive.disorder + intraparenchymal, 
                data = chart.min)

sens_min = sensemakr(model = lin_min, treatment = "intervention.group",
                     benchmark_covariates = c("age", "gcs"), 
                     kd = 1:3)

# Display Statistics of Sensitivity Analysis

ovb_minimal_reporting(sens_min, format = "latex")


# label.bump.x moves the points, change accordingly
plot(sens_min, xlab = "Partial R² of Confounders with Minimally Invasive Surgery",
     ylab = "Partial R² of Confounders with Mortality", cex.lab = 1.25, 
     family = "Arial", label.bump.x = 1)

# Craniotomy and Craniectomy Sensitivity Analysis
# _____________________________________________________________________________

# Linear regression model
lin_cranio = lm(died.before.90.days ~ intervention.group + age + 
               gcs + wbc + platelet + hematocrit + glucose + 
               creatinine + sodium + potassium + diabetes + anemia + 
               liver.disease + depressive.disorder + intraparenchymal, 
             data = chart.cranio)

sens_cranio = sensemakr(model = lin_cranio, treatment = "intervention.group",
                     benchmark_covariates = c("age"), 
                                              kd = 1:3)

# Display Statistics of Sensitivity Analysis
ovb_minimal_reporting(sens_cranio, format = "latex")

# label.bump.x moves the points, change accordingly
plot(sens_cranio, xlab = "Partial R² of Confounders with Craniotomy/Craniectomy",
     ylab = "Partial R² of Confounders with Mortality", cex.lab = 1.25, 
     family = "Arial", label.bump.x = 1)
