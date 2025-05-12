################################################################################
# Title: Internal Consistency
################################################################################

library(tidyverse)
library(here)
library(labelled)
library(skimr)
library(geomtextpath)
library(effectsize)
library(vcdExtra)
library(ggpubr)
library(viridis)
library(gt)
library(psych)
library(officer)
library(flextable)


# Load new data 
mydata <- readRDS(here("data_ori", "name_data.rds"))

# Remove label

nolabel_data <- remove_val_labels(mydata) |> 
  mutate(Edad_Recodificada = case_when(
    Edad_Recodificada == 15 ~ 2,
    TRUE ~ Edad_Recodificada))


# Select only the data for worry about risk perceptions
data <- nolabel_data |> 
  select(2:7, 14, 40:57, 59:80, 83:103, 106:123)

# Select variables to perform alpha
dataset_list<- list(
  
worry_cyber <- data[, c(8:16)],

risk_cyber <- data[, c(17:25)],

worry_bully <- data[, c(26:36)],

risk_bully <- data[, c(37:47)],

worry_partner <- data[, c(51:59)],

risk_partner <- data[, c(60:68)],

worry_sex_vio <- data[, c(69:77)],

risk_sex_vio <- data[, c(78:86)])


# Perform alpha for each scale with bootstrapped 95 % CI with 1000 resamples

internal_consistency <- lapply(dataset_list,  alpha, n.iter = 1000)

name_scales <- list("worry_cyber", "risk_cyber", "worry_bully" , "risk_bully",
                    "worry_partner", "risk_partner", "worry_sex_vio", "risk_sex_vio")

internal_cons_resul_list <- list()

for (i in seq(length(internal_consistency))) {
  
  
  internal_cons_resul_list[[i]] <- internal_consistency[[i]]$total |> 
    mutate(lower_ci = internal_consistency[[i]]$boot.ci[[1]],
           upper_ci = internal_consistency[[i]]$boot.ci[[3]],
           name_scale = name_scales[[i]]) |> 
    select(name_scale, raw_alpha, std.alpha, lower_ci, upper_ci)
    
  
}


internal_consis_resul_tb <- bind_rows(internal_cons_resul_list)

# Print table in docx

if (!file.exists(here("outcomes", "internal_consist_result.docx"))) {
  
  # Convert the data frame to a flextabe
  internal_consis_fxtb <- flextable(internal_consis_resul_tb)
  
  # create a word document and add the table
  internal_consis_fxtb <- read_docx() |> 
    body_add_flextable(value = internal_consis_fxtb) 
  
  # Save the table 
  print(internal_consis_fxtb, target = here("outcomes", "internal_consist_result.docx"))
  
} else {
  print("Table already exists!")
  
}

