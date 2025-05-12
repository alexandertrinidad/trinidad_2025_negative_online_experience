################################################################################
# Title: Exploring missing data

# Loading packages and data -----------------------------------------------

library(tidyverse)
library(here)
library(sjPlot)
library(labelled)
library(VIM)
library(naniar)


# Load new data.

mydata <- readRDS(here("data_ori", "name_data.rds"))

# Remove label

nolabel_data <- remove_val_labels(mydata) 


# Checking duplicated cases ------------------------------------------------
# Remove character variables to check duplicated rows with numeric values
only_num <- nolabel_data |> 
  select(where(is.numeric))

glimpse(only_num)

# Check if completley cases

nrow(distinct(only_num))

# Identify same ID
dupli_id <- nolabel_data |> 
  filter(duplicated(ID)) |> 
  select(ID)

# Check if same ID 
same_id <- nolabel_data |> 
  filter(ID %in% dupli_id$ID)

# After inspection it seems to be a coding error. Rows are completely different.
# Check if duplicates in the CODIGO column

table(duplicated(nolabel_data$CODIGO)) # No duplicated. This is the unique key


# Missing Data Patterns ---------------------------------------------------


# Completeness test
table(stats::complete.cases(only_num)) # onlny 24 cases are completed


# Missing Data Patterns ---------------------------------------------------

# Demographic variables
demogr_data <- nolabel_data[ , c(3,4, 5:7)]
summary(aggr(demogr_data, sortVar=TRUE))$combinations
matrixplot(demogr_data)

# Cyber-crime 
## Worry:
freq_worry_cyber <- nolabel_data[ , c(40:48)]
summary(aggr(freq_worry_cyber, sortVar=TRUE))$combinations


## Risk perception
percb_risk_cyber <- nolabel_data[ , c(49:57)]
summary(aggr(percb_risk_cyber, sortVar=TRUE))$combinations

# Cyberbullying
## Worry:
freq_worry_bully <- nolabel_data[ , c(59:69)]
summary(aggr(freq_worry_bully, sortVar=TRUE))$combinations

## Risk perception
percb_risk_bully <- nolabel_data[, c(70:80)]
summary(aggr(percb_risk_bully, sortVar=TRUE))$combinations

# Dating violence
## Worry
freq_worry_partner <- nolabel_data[ , c(86:94)]
summary(aggr(freq_worry_partner, sortVar=TRUE))$combinations

## Risk perception
percb_risk_partner <- nolabel_data[ , c(95:103)]
summary(aggr(percb_risk_partner, sortVar=TRUE))$combinations

# Sexual violence
## Worry
freq_worry_sex_vio <- nolabel_data[ , c(106:114)]
summary(aggr(freq_worry_sex_vio, sortVar=TRUE))$combinations

## Risk perception
percb_worry_sex_vio <- nolabel_data[ , c(115:123)]
summary(aggr(percb_worry_sex_vio, sortVar=TRUE))$combinations


# Missing Data Evaluation -------------------------------------------------
# Little's MCAR test

# Selecting only the variables of interest
data_study <- nolabel_data |> 
  select(c(2,3,4, 5:7, 40:57, 59:80, 86:103, 106:123))

mcar_test(data = data_study) 

