################################################################################
# Title: Prevalence of online negative experience and perceived risks.
################################################################################
# Load packages
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
library(officer)
library(cowplot)

# Load new data 
mydata <- readRDS(here("data_ori", "name_data.rds"))

# Remove label

nolabel_data <- remove_val_labels(mydata) |> 
  mutate(Edad_Recodificada = case_when(
    Edad_Recodificada == 15 ~ 2,
    TRUE ~ Edad_Recodificada))


# Data preparation --------------------------------------------------------


# Select only the data for worry about risk perceptions
data <- nolabel_data |> 
  select(2:7, 14, 40:57, 59:80, 83:103, 106:123)

# Create scale of Worry and risk perception
# by averaging the values (to deal with missing values)
data_scales <- data |> 
  rowwise() |> 
  mutate(freq_worry_cyber = mean(c_across(c(8:16)), na.rm = TRUE),
         percv_risk_cyber = mean(c_across(c(17:25)), na.rm = TRUE),
         freq_worry_bully = mean(c_across(c(26:36)), na.rm = TRUE),
         percv_risk_bully = mean(c_across(c(37:47)), na.rm = TRUE),
         freq_worry_partner = mean(c_across(c(51:59)), na.rm = TRUE),
         percv_risk_partner = mean(c_across(c(60:68)), na.rm = TRUE),
         freq_worry_sex_vio = mean(c_across(c(69:77)), na.rm = TRUE),
         percv_worry_sex_vio = mean(c_across(c(78:86)), na.rm = TRUE)) |> 
  ungroup() |> 
  select(1:7, 87:94) |> 
  mutate(Id.sexual = as.factor(Id.sexual),
         Edad_Recodificada = as.ordered(Edad_Recodificada),
         Curso_escolar = as.ordered(Curso_escolar),
         Centro_escolar = as.factor(Centro_escolar),
         sit_negat = case_when(
           A.Sit.negativa == 2 ~ 0,
           TRUE ~ A.Sit.negativa),
         age_recoded = case_when(
           Edad == 12 ~ "12-14",
           Edad == 13 ~ "12-14",
           Edad == 14 ~ "12-14",
           Edad == 15 ~ "15-18",
           Edad == 16 ~ "15-18",
           Edad == 17 ~ "15-18",
           Edad == 18 ~ "15-18"))
  
# Descriptive
# skim(data_scales) # Activate to explore the descriptive

# Re-coding high-perceived worry and risk (mean +/- sd)
# Perceived risk 3 or more is perceived risk or high risk, meanwhile in frequency
# we use the mean +/- sd to determine high frequency
data_scales <- data_scales |> 
  mutate(high_freq_worry_cyber = case_when(
    freq_worry_cyber >= (mean(data_scales$freq_worry_cyber, na.rm = TRUE) +
                           sd(data_scales$freq_worry_cyber, na.rm = TRUE))  ~ 1,
    freq_worry_cyber < (mean(data_scales$freq_worry_cyber, na.rm = TRUE) +
                          sd(data_scales$freq_worry_cyber, na.rm = TRUE)) ~ 0),
    di_percv_risk_cyber = case_when(
      percv_risk_cyber >= 3 ~ 1,
      percv_risk_cyber < 3 ~ 0),
    high_freq_worry_bully = case_when(
      freq_worry_bully >= (mean(data_scales$freq_worry_bully, na.rm = TRUE) +
                             sd(data_scales$freq_worry_bully, na.rm = TRUE))  ~ 1,
      freq_worry_bully < (mean(data_scales$freq_worry_bully, na.rm = TRUE) +
                            sd(data_scales$freq_worry_bully, na.rm = TRUE)) ~ 0),
    di_percv_risk_bully = case_when(
      percv_risk_bully >= 3 ~ 1,
      percv_risk_bully < 3 ~ 0),
    high_freq_worry_partner = case_when(
      freq_worry_partner >= (mean(data_scales$freq_worry_partner, na.rm = TRUE) +
                             sd(data_scales$freq_worry_partner, na.rm = TRUE))  ~ 1,
      freq_worry_partner < (mean(data_scales$freq_worry_partner, na.rm = TRUE) +
                            sd(data_scales$freq_worry_partner, na.rm = TRUE)) ~ 0),
    di_percv_risk_partner = case_when(
      percv_risk_partner >= 3 ~ 1,
      percv_risk_partner < 3 ~ 0),
    high_freq_worry_sex_vio = case_when(
      freq_worry_sex_vio >= (mean(data_scales$freq_worry_sex_vio, na.rm = TRUE) +
                               sd(data_scales$freq_worry_sex_vio, na.rm = TRUE))  ~ 1,
      freq_worry_sex_vio < (mean(data_scales$freq_worry_sex_vio, na.rm = TRUE) +
                              sd(data_scales$freq_worry_sex_vio, na.rm = TRUE)) ~ 0),
    di_percv_risk_sex_vio = case_when(
      percv_worry_sex_vio >= 3 ~ 1,
      percv_worry_sex_vio < 3 ~ 0))


# Prevalence of Negative Experience ---------------------------------------

# Proportions chi-squared test against 0.5 probability 
prevalence_prop <- prop.test(table(data_scales$sit_negat)[[2]], 
                             table(is.na(data_scales$sit_negat))[[1]],
                             correct = FALSE)

# Prevalence 
prevalence_value <- prevalence_prop$estimate[[1]]

# Cohen's h' with bootstrapping 
cohens_h_ci_bootstrap <- function(data, var1, p1, p0, n_bootstrap = 1000,
                                  confidence_level = 0.95, plot = FALSE) {
  # Function to compute Cohen's h
  compute_h_prime <- function(p1, p0) {
    h_ori <- abs(2 * asin(sqrt(p1)) - 2 * asin(sqrt(p0)))
    
    h_ori * sqrt(2)
  }
  
  # Calculate original Cohen's h
  h_prime <- compute_h_prime(p1, p0)
  
  var1 <- ensym(var1)
  
  set.seed(123)
  
  if (plot) {
    
    # Bootstrap resampling
    bootstrap_hs <- replicate(n_bootstrap, {
      sample_data <- slice_sample(data, n = nrow(data), replace = TRUE)
      prevalence_prop <- prop.test(table(sample_data[[var1]])[[2]], table(is.na(sample_data[[var1]]))[[1]], correct = FALSE)
      sample_p1 <- prevalence_prop$estimate[[1]]
      compute_h_prime(sample_p1, p0)
    })
    
    dt_hs <- as.data.frame(bootstrap_hs)
    
    # Calculate confidence intervals
    alpha <- 1 - confidence_level
    ci_lower <- quantile(bootstrap_hs, probs = alpha / 2, na.rm = TRUE)
    ci_upper <- quantile(bootstrap_hs, probs = 1 - alpha / 2, na.rm = TRUE)
    
    # Return results
    myhbootlist <- list(h_prime = h_prime, ci_lower = ci_lower, ci_upper = ci_upper)
    
    myplot <- ggplot(data = dt_hs, aes(x = bootstrap_hs)) +
      geom_density(fill = "blue", alpha = 0.3) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      labs(x = paste("Bootstrapped Cohen's h' with ", n_bootstrap, " resamples"),
           y = "Density") +
      geom_vline(xintercept = myhbootlist$h_prime, color = "darkred") +
      geom_vline(xintercept = myhbootlist$ci_lower, color = "darkgreen", linetype = "dashed") +
      geom_vline(xintercept = myhbootlist$ci_upper, color = "darkgreen", linetype = "dashed") +
      annotate("text", x = myhbootlist$h_prime, y = Inf, label = "Cohen's h'", color = "darkred", vjust = -1.5) +
      annotate("text", x = myhbootlist$ci_lower, y = Inf, label = "Lower 95% CI", color = "darkgreen", vjust = -1.5, angle = 90, hjust = 1.2) +
      annotate("text", x = myhbootlist$ci_upper, y = Inf, label = "Upper 95% CI", color = "darkgreen", vjust = -1.5, angle = 90, hjust = 1.2)
    
    return(myplot)
    
  } else {
    
    # Bootstrap resampling
    bootstrap_hs <- replicate(n_bootstrap, {
      sample_data <- slice_sample(data, n = nrow(data), replace = TRUE)
      prevalence_prop <- prop.test(table(sample_data[[var1]])[[2]], table(is.na(sample_data[[var1]]))[[1]], correct = FALSE)
      sample_p1 <- prevalence_prop$estimate[[1]]
      compute_h_prime(sample_p1, p0)
    })
    
    # Calculate confidence intervals
    alpha <- 1 - confidence_level
    ci_lower <- quantile(bootstrap_hs, probs = alpha / 2, na.rm = TRUE)
    ci_upper <- quantile(bootstrap_hs, probs = 1 - alpha / 2, na.rm = TRUE)
    
    # Return results
    list(h_prime = h_prime, ci_lower = ci_lower, ci_upper = ci_upper)
    
  }
}


my_h_plot <- cohens_h_ci_bootstrap(data_scales, var1 = sit_negat, prevalence_value, 0.5, plot = TRUE)

# Print Figure 1
ggsave("TrinidadFigure1.tiff", dpi = 300)
dev.off()


# Cohen's h' with bootstrapped 95 % CI 
my_h_ci <- cohens_h_ci_bootstrap(data_scales, var1 = sit_negat,
                                p1 = prevalence_value, 
                                p0 = 0.5,
                                plot = FALSE)

# Prevalence of Negative Experience According to Demographics -------------

# Recoding variables
data_scales <- data_scales |> 
  mutate(id_sex_labs = case_when(
    Id.sexual == 1 ~ "Female",
    Id.sexual == 2 ~ "Male",
    TRUE ~ NA), # Remove non cisgender because there are only 15,
    cat_h_freq_worry_cyber = case_when(
      high_freq_worry_cyber == 1 ~ "Yes", 
      high_freq_worry_cyber == 0 ~ "No", 
      TRUE ~ NA),
    cat_percv_risk_cyber = case_when(
      di_percv_risk_cyber == 1 ~ "Yes",
      di_percv_risk_cyber == 0 ~ "No",
      TRUE ~ NA),
    cat_h_freq_worry_bully = case_when(
      high_freq_worry_bully == 1 ~ "Yes",
      high_freq_worry_bully == 0 ~ "No",
      TRUE ~ NA),
    cat_percv_risk_bully = case_when(
      di_percv_risk_bully == 1 ~ "Yes",
      di_percv_risk_bully == 0 ~ "No",
      TRUE ~ NA),
    cat_h_freq_worry_partner = case_when(
      high_freq_worry_partner == 1 ~ "Yes",
      high_freq_worry_partner == 0 ~ "No",
      TRUE ~ NA),
    cat_percv_risk_partner = case_when(
      di_percv_risk_partner == 1 ~ "Yes",
      di_percv_risk_partner == 0 ~ "No",
      TRUE ~ NA),
    cat_h_freq_worry_sex_vio = case_when(
      high_freq_worry_sex_vio == 1 ~ "Yes",
      high_freq_worry_sex_vio == 0 ~ "No",
      TRUE ~ NA),
    cat_percv_risk_sex_vio = case_when(
      di_percv_risk_sex_vio == 1 ~ "Yes",
      di_percv_risk_sex_vio == 0 ~ "No",
      TRUE ~ NA))

# Remove non-binary
data_scale_no_non_binary <- data_scales |> 
  dplyr::filter(!is.na(id_sex_labs))

# Prevalence of Negative Online Experience According to Demographic Variables
mylist <- vector(mode = "list", length = 2)
index <- 1

for (i in c(17,26)) {
  
  myunique_level <- unique(na.omit(data_scale_no_non_binary[i]))[[1]]
  mycontig_table <- table(data_scale_no_non_binary[[i]], data_scale_no_non_binary[[16]])
  mypercentages <- prop.table(mycontig_table, margin = 1)
  chi_square_result <- chisq.test(mycontig_table)
  
if(length(unique(na.omit(data_scale_no_non_binary[[i]]))) > 2) {
  
  mycramerv <- cramers_v(mycontig_table, adjust = FALSE)
  
  mytable <- tibble("Varibles" = myunique_level, 
                    "Negative Experience (%)" = c(mypercentages[1,2]*100, mypercentages[2,2]*100,  mypercentages[3,2]*100),
                    "\u03C7\u00B2" = c(chi_square_result$statistic[[1]], NA, NA),
                    "p-value" = c(chi_square_result$p.value, NA, NA),
                    "Phi / V" = c(round(mycramerv$Cramers_v, 3), NA, NA),
                    "95 % CI" = c(paste("[",round(mycramerv$CI_low, 2), ",", round(mycramerv$CI_high, 2), "]"), NA, NA))
} else {
  
  myphi <- effectsize::phi(mycontig_table, adjust = FALSE, alternative = "two.sided")
  
  mytable <- tibble("Variables" = myunique_level, 
                    "Negative Experience (%)" = c(mypercentages[1,2]*100, mypercentages[2,2]*100),
                    "\u03C7\u00B2" = c(chi_square_result$statistic[[1]], NA),
                    "p-value" = c(chi_square_result$p.value, NA),
                    "\u03C6" = c(round(myphi$phi, 3), NA),
                    "95 % CI" = c(paste("[",round(myphi$CI_low, 2), ",", round(myphi$CI_high, 2), "]"), NA))
  
}
 
  
  
  
  mylist[[index]] <- mytable
  index <- index + 1
  
}

# All results in a table 
combined_results <- bind_rows(mylist)

if (!file.exists(here("outcomes", "table_1_prev_demogra.docx"))) {

# Convert the data frame to a flextabe
ft_prev_demogra <- flextable(combined_results)

# create a word document and add the table
doc_prev_demogra <- read_docx() |> 
  body_add_flextable(value = ft_prev_demogra) 

# Save the table 
print(doc_prev_demogra, target = here("outcomes", "table_1_prev_demogra.docx"))

} else {
  print("Table already exists!")
  
}

# Risk Perception and High Frequency Worry per sexual Id.----------------

  
mylist_high_perc <- vector(mode = "list", length = 8)
index <- 1

for (i in 27:34) {
  
  myvariable_name <- names(data_scale_no_non_binary[i])
  
  myvariab_percen <- (sum(data_scale_no_non_binary[i] == "Yes", na.rm = TRUE)/
                        sum(complete.cases(data_scale_no_non_binary[i]))) * 100
  
  mycontig_table <- table(data_scale_no_non_binary[[i]], data_scale_no_non_binary[[26]])
  mypercentages <- prop.table(mycontig_table, margin = 2)
  chi_square_result <- chisq.test(mycontig_table)
  
 
    
    myphi <- effectsize::phi(mycontig_table, adjust = FALSE, alternative = "two.sided")
    
    mytable <- tibble("High..." = myvariable_name, 
                      "N" = sum(complete.cases(data_scale_no_non_binary[i])),
                      "Total (%)" = round(myvariab_percen, 2),
                      "Female %" =  round(mypercentages[2,1]*100, 2),
                      "Male %" = round(mypercentages[2,2]*100, 2),
                      "\u03C7\u00B2" = chi_square_result$statistic[[1]],
                      "p-value" = chi_square_result$p.value,
                      "\u03C6" = round(myphi$phi, 3),
                      "95 % CI" = paste("[",round(myphi$CI_low, 2), ",", round(myphi$CI_high, 2), "]"))

  mylist_high_perc[[index]] <- mytable
  index <- index + 1
  
}

# All results in a table 1A
combined_results_perc <- bind_rows(mylist_high_perc) 
if (!file.exists(here("outcomes", "table_1A_percv_worry_gender.docx"))) {
  
  # Convert the data frame to a flextabe
  wor_percv_gender <- flextable(combined_results_perc)
  
  # create a word document and add the table
  doc_wor_percv_gender <- read_docx() |> 
    body_add_flextable(value = wor_percv_gender) 
  
  # Save the table 
  print(doc_wor_percv_gender, target = here("outcomes", "table_1A_percv_worry_gender.docx"))
  
} else {
  print("Table already exists!")
  
}

# Divide the frequency of worry and risk in two tables
h_freq_worry <- combined_results_perc |> 
  filter(startsWith(combined_results_perc[[1]], "cat_h_freq"))

h_freq_worry_long <- pivot_longer(h_freq_worry, cols = c(4, 5), 
                                  names_to = "id_sex", 
                                  values_to = "percentage")

h_freq_worry_long_selct <- h_freq_worry_long |> 
  select(1, 5, 6, 8, 9) |> 
  dplyr::rename(type_of_off = names(h_freq_worry_long)[1]) |> 
  mutate(id_sex = case_when(
    id_sex == "Female %" ~ "Female",
    id_sex == "Male %" ~ "Male")) 


# Plot 1
plot1 <- ggplot(h_freq_worry_long_selct, aes(x = type_of_off, 
                                             y = percentage,
                                             fill = id_sex)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_signif(
    mapping = aes(x = type_of_off),
    y_position = c(26, 23, 21, 23),
    xmin = c(0.5, 1.5, 2.5, 3.5),
    xmax = c(1.4, 2.4, 3.4, 4.4),
    annotation = c(paste0("\u03C6" , " = ", h_freq_worry_long_selct[[3,3]] ," ***"),
                   paste0("\u03C6" , " = ", h_freq_worry_long_selct[[1,3]] ," ***"), 
                   paste0("\u03C6" , " = ", h_freq_worry_long_selct[[5,3]] ," ***"), 
                   paste0("\u03C6" , " = ", h_freq_worry_long_selct[[7,3]] ," ***"))) +
  geom_text(
    aes(label = format(round(percentage, 2), nsmall = 1)),
    position = position_dodge(width = 0.9),
    colour = "black",
    hjust = 0.65,
    vjust = 4) +
  scale_fill_viridis(discrete = TRUE,
                     option = "D",
                     begin = 0.80,
                     end = 1,
                     alpha = 1) +
  scale_x_discrete(
    labels = c("cat_h_freq_worry_bully" = "Bullying",
               "cat_h_freq_worry_cyber" = "Cyber", 
               "cat_h_freq_worry_partner"= "Dating Violence",
               "cat_h_freq_worry_sex_vio" = "Sexual Violence")) +
  labs(title = "High Frequency of Worry about Online Victimization",
       x = "",
       y = "% of High Frequency of Worry",
       fill = "",
       cex.main = 2) +
  theme_minimal()  +
  theme(legend.position = "bottom")
  


# Select only the risk perception
risk_percep <- combined_results_perc |> 
  filter(startsWith(combined_results_perc[[1]], "cat_percv_ri"))

risk_percep_long <- pivot_longer(risk_percep, cols = c(4, 5), 
                                  names_to = "id_sex", 
                                  values_to = "percentage")

risk_percep_long_selct <- risk_percep_long |> 
  select(1, 5, 6, 8, 9) |> 
  dplyr::rename(type_of_off = names(risk_percep_long)[1]) |> 
  mutate(id_sex = case_when(
    id_sex == "Female %" ~ "Female",
    id_sex == "Male %" ~ "Male")) 
  
# Plot 2 
plot2 <- ggplot(risk_percep_long_selct, aes(x = type_of_off, 
                                                   y = percentage,
                                                   fill = id_sex)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_signif(
    mapping = aes(x = type_of_off),
    y_position = c(11, 6, 5, 9),
    xmin = c(0.5, 1.5, 2.5, 3.5),
    xmax = c(1.4, 2.4, 3.4, 4.4),
    annotation = c(paste0("\u03C6" , " = ", risk_percep_long_selct[[3,3]] ," *"), 
                   paste0("\u03C6" , " = ", risk_percep_long_selct[[1,3]] ," ***"),
                   paste0("\u03C6" , " = ", risk_percep_long_selct[[5,3]] ,""), 
                   paste0("\u03C6" , " = ", risk_percep_long_selct[[7,3]] ," **"))) +
  geom_text(
    aes(label = format(round(percentage, 2), nsmall = 1)),
    position = position_dodge(width = 0.9),
    colour = "black",
    hjust = 0.65,
    vjust = 1) +
  scale_fill_viridis(discrete = TRUE,
                     option = "D",
                     begin = 0.60,
                     end = 0.3,
                     alpha = 1) +
  scale_y_continuous(limits = c(0, 25)) +
  scale_x_discrete(
    labels = c("cat_percv_risk_bully" = "Bullying", 
               "cat_percv_risk_cyber" = "Cyber",
               "cat_percv_risk_partner" = "Dating Violence", 
               "cat_percv_risk_sex_vio" = "Sexual Violence")) +
  labs(title = "Risk Perception of Online Victimization",
       x = "",
       y = "% of Risk Perception",
       fill = "",
       cex.main = 2) +
  theme_minimal()  +
  theme(legend.position = "bottom")

# Figure 
plot_grid(plot1, plot2, labels = "AUTO")

ggsave(here("submission","TrinidadFigure2.tiff"), scale = 2 ,limitsize = FALSE, dpi = 300)
dev.off()



# Three way interactions Frequency Worry
high_worry_tabs_list <- vector(mode = "list", 4)
index <- 1

for (i in c(27,29,31,33)) {
  
  high_worry_tabs <- xtabs(~ data_scale_no_non_binary[[i]] + data_scale_no_non_binary[[26]] +
                             data_scale_no_non_binary[[16]], #  + data_scale_no_non_binary[[17]], 
                           data = data_scale_no_non_binary)
  
  
  high_worry_tabs_list[[index]] <- high_worry_tabs
  index <- index + 1
}


contingency_tabs <- list()

for (i in seq_along(high_worry_tabs_list)) {
  
  # Log-linear model per type of crime 
  type_of_offense <- as.data.frame(high_worry_tabs_list[[i]]) 
  
  variable_position <- c(27,29,31,33)
  
  varname <- names(data_scale_no_non_binary)[variable_position[i]]
    
    # Old and new variable names
    old_names <- names(type_of_offense)
    new_names <- c(varname, "id_sex", "negat_exp", "Freq")
    
    # Create a named vector for renaming
    rename_vector <- setNames(new_names, old_names)
    
    # Rename columns using rename_with and a function
    renamed_variables <- type_of_offense  |> 
      rename_with(~ rename_vector[.x], all_of(old_names))
    
    
    contingency_tabs[[paste0("contingency_tab_", i)]] <- renamed_variables
    
  }


# Assign each data frame in the list to the environment


# Log linear models Worry -------------------------------------------------------
best_model_typ_worry <- list()

for (i in seq_along(contingency_tabs)) {
  

dataset_tabs <- contingency_tabs[[i]]

varname <- dataset_tabs[[1]]

# Select variable of reference
dataset_tabs$id_sex = relevel(dataset_tabs$id_sex, ref = "Female")
dataset_tabs$negat_exp = relevel(dataset_tabs$negat_exp, ref = "1")


# Fit satuarted model (XYZ)
sat.model <- glm(Freq ~ varname * id_sex * negat_exp, 
                     family = poisson(), data = dataset_tabs)

# Fit complete independence (X, Y, Z)

compl.indep <- glm(Freq ~ varname + id_sex + negat_exp, 
                   family = poisson(), data = dataset_tabs)

# Joint independence (X, YZ)
joint.indep <- glm(Freq ~ varname + id_sex + negat_exp + 
                     negat_exp*varname, 
                   family = poisson(), data = dataset_tabs)

# Conditional independence (XZ, YZ)
condi.indep <- glm(Freq ~ varname + id_sex + negat_exp + 
                     negat_exp*varname + negat_exp*id_sex, 
                   family = poisson(), data = dataset_tabs)

# Homogeneous Association (XZ, YZ, XY)

homog.associ <- glm(Freq ~ (varname + id_sex + negat_exp)^2, 
                    family = poisson(), data = dataset_tabs)

# Model list

models_list <- list(sat.model, compl.indep, joint.indep, condi.indep, homog.associ)


# Model fit

models_fit <- vcdExtra::LRstats(glmlist(sat.model,
                                        compl.indep,
                                        joint.indep, 
                                        condi.indep,
                                        homog.associ))

models_fit_dt <- as.data.frame(models_fit)

if(!file.exists(here("outcomes", paste0("best_model_fit_", i, ".docx")))) {
  
  models_fit_fx <- flextable(models_fit_dt)
  
  # create a word document and add the table
  doc_models_fit_fx <- read_docx() |> 
    body_add_flextable(value = models_fit_fx) 
  
  # Save the table 
  print(doc_models_fit_fx, target = here("outcomes", paste0("best_model_fit_", i, ".docx")))
  
  
}



# Select the model based on AIC

lowest_AIC_index <- which.min(models_fit$AIC)

best_model <- models_list[[lowest_AIC_index]]

best_model_typ_worry[[paste0("best_model_", i)]] <- best_model

}

# Create a data frame with the two-ways or three ways interactions

result_transfor_tb <- function(index, contingency_tb, best_model) {
  
  # Name of the type of victimization 
  type_of_victimization <- names(contingency_tb[[index]])[1]
  
  #three two-ways interacitions
  model_tidy <- best_model[[index]] |> broom::tidy() #
  
  myinteraction_transformed <- list()
  
  for (i in 5:nrow(model_tidy)) {
    
    
    interaciton_type <- model_tidy$term[i] 
    coef_interaction <- model_tidy$estimate[i] 
    se_interaction <- model_tidy$std.error[i]
    
    # Transforming standard error of cohen's d from logOrs 
    se_logodds_to_cohens_d <- function(se_lor) {
      transformation_factor <- sqrt(3)  / pi
      se_d <- se_lor * transformation_factor
      return(se_d)
    }
    
    se_cohens_d <-se_logodds_to_cohens_d(se_interaction)
    
    # Transfomr logOdds to Cohens'd
    cohens_d <- effectsize::logoddsratio_to_d(coef_interaction)
    
    # Calculating CIs
    low_ci_cohens_d <- function(d, se, z){
      
      ci_low <- d - z*se
      return(ci_low)
    }
    
    ci_low <- low_ci_cohens_d(cohens_d, se_cohens_d, 1.96)
    
    high_ci_cohens_d <- function(d, se, z){
      
      ci_high <- d + z*se
      return(ci_high)
    }
    
    ci_high <- high_ci_cohens_d(cohens_d, se_cohens_d, 1.96)
    
    
    # Create a tibble with the data
    transformed_results <- tibble(
      type_of_victimization = type_of_victimization,
      models_interaction = interaciton_type,
      cohens_d = cohens_d,
      ci_low = ci_low,
      ci_high = ci_high)
    
    myinteraction_transformed[[i]] <- transformed_results
    
  }
  
  transformed_interaction_tb <- bind_rows(myinteraction_transformed)
}

transfor_result_worry <- list()

for (i in seq_along(best_model_typ_worry)) {
  
  transfor_result_worry[[i]] <- result_transfor_tb(i, contingency_tabs, best_model_typ_worry)
 }

transfor_result_worry_all <- bind_rows(transfor_result_worry)

# Check the unique values in type_of_victimization
actual_levels <- unique(transfor_result_worry_all$type_of_victimization)

# Changing order 
transfor_result_worry_all$type_of_victimization = factor(transfor_result_worry_all$type_of_victimization, 
                                                          levels = actual_levels)
# Pattern to remove
pattern_remove <- "cat_h_freq_worry_"
# Rename interactions 
transfor_result_worry_all <- transfor_result_worry_all |> 
  mutate(worry_offense = gsub(pattern_remove, "", type_of_victimization),
         interactions = case_when(
           worry_offense == "cyber" &
             models_interaction == "varnameYes:id_sexMale" ~ "Worry about Cyber: High worry x Male",
           worry_offense == "cyber" &
             models_interaction == "varnameYes:negat_exp0" ~ "Worry about Cyber: High worry x No Negative Exp.",
           worry_offense == "cyber" &
             models_interaction == "id_sexMale:negat_exp0" ~ "Worry about Cyber: Male x No Negative Exp.",
           worry_offense =="cyber" &
             models_interaction == "varnameYes:id_sexMale:negat_exp0" ~ "Worry about Cyber: High worry x Male x No Negative Exp.",
           worry_offense == "bully" &
             models_interaction == "varnameYes:id_sexMale" ~ "Worry about Bullying: High worry x Male",
           worry_offense == "bully" &
             models_interaction == "varnameYes:negat_exp0" ~ "Worry about Bullying: High worry x No Negative Exp.",
           worry_offense == "bully" &
             models_interaction == "id_sexMale:negat_exp0" ~ "Worry about Bullying: Male x No Negative Exp.",
           worry_offense == "partner" &
             models_interaction == "varnameYes:id_sexMale" ~ "Worry about Dating Violence: High worry x Male",
           worry_offense == "partner" &
             models_interaction == "varnameYes:negat_exp0" ~ "Worry about Dating Violence: High worry x No Negative Exp.",
           worry_offense == "partner" &
             models_interaction == "id_sexMale:negat_exp0" ~ "Worry about Dating Violence: Male x No Negative Exp.",
           worry_offense == "sex_vio" &
             models_interaction == "varnameYes:id_sexMale" ~ "Worry about Sexual Harassment: High worry x Male",
           worry_offense == "sex_vio" &
             models_interaction == "varnameYes:negat_exp0" ~ "Worry about Sexual Harassment: High worry x No Negative Exp.",
           worry_offense == "sex_vio" &
             models_interaction == "id_sexMale:negat_exp0" ~ "Worry about Sexual Harassment: Male x No Negative Exp."))

# Save the table as supplementary Table S3
if(!file.exists(here("outcomes", paste0("log_lin_res_worry.docx")))) {
  
  transfor_result_worry_fx <- flextable(transfor_result_worry_all)
  
  # create a word document and add the table
  doc_transfor_result_worry <- read_docx() |> 
    body_add_flextable(value = transfor_result_worry_fx) 
  
  # Save the table 
  print(doc_transfor_result_worry, target = here("outcomes", paste0("log_lin_res_worry.docx")))
  
  
}

# Define a function to extract the part of the string before the underscore
extract_after_underscore <- function(label) {
  sub(".*_", "", label)
}

# Define a function to extract the part of the string before the underscore
extract_after_colon <- function(label) {
  sub(".*:", "", label)
}

# Define a function to extract the part of the string before the underscore
extract_before_colon <- function(label) {
  sub(":.*", "", label)
}

transfor_result_worry_all <- transfor_result_worry_all |> 
  mutate(interactions_gen = extract_after_colon(interactions),
         high_worry = extract_before_colon(interactions))

# Forestplot worry about points = interactions y =  models. . 

forest_plot1 <- transfor_result_worry_all |> 
  ggplot(aes(y = high_worry, x = cohens_d, color = interactions_gen, fill = interactions_gen)) +
  geom_point(shape = 21, size = 4, position = position_dodge(width = 0.75)) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), 
                 position = position_dodge(width = 0.75), linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "", 
       x = "Cohen's d",
       fill = "",
       color = "") +
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        text = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 2), 
         color = guide_legend(nrow = 2)) +
  scale_fill_viridis(discrete = TRUE,
                     option = "D",
                     begin = 0,
                     end = 0.9) +
  scale_color_viridis(discrete = TRUE,
                      option = "D",
                      begin = 0,
                      end = 0.9) 

# Three way interactions Risk
risk_tabs_list <- vector(mode = "list", 4)
index <- 1

for (i in c(28,30,32,34)) {
  
  high_risk_tabs <- xtabs(~ data_scale_no_non_binary[[i]] + data_scale_no_non_binary[[26]] +
                             data_scale_no_non_binary[[16]], 
                           data = data_scale_no_non_binary)
  
  
  risk_tabs_list[[index]] <- high_risk_tabs
  index <- index + 1
}


contingency_tabs_risk <- list()

for (i in seq_along(risk_tabs_list)) {
  
  # Log-linear model per type of crime 
  type_of_offense <- as.data.frame(risk_tabs_list[[i]]) 
  
  variable_position <- c(28,30,32,34)
  
  varname <- names(data_scale_no_non_binary)[variable_position[i]]
  
  # Old and new variable names
  old_names <- names(type_of_offense)
  new_names <- c(varname, "id_sex", "negat_exp", "Freq")
  
  # Create a named vector for renaming
  rename_vector <- setNames(new_names, old_names)
  
  # Rename columns using rename_with and a function
  renamed_variables <- type_of_offense  |> 
    rename_with(~ rename_vector[.x], all_of(old_names))
  
  
  contingency_tabs_risk[[paste0("contingency_tab_risk", i)]] <- renamed_variables
  
}

# Log linear models Risk perception -------------------------------------------------------
best_model_typ_risk <- list()

for (i in seq_along(contingency_tabs_risk)) {
  
  
  dataset_tabs <- contingency_tabs_risk[[i]]
  
  varname <- dataset_tabs[[1]]
  
  # Select variable of reference
  dataset_tabs$id_sex = relevel(dataset_tabs$id_sex, ref = "Female")
  dataset_tabs$negat_exp = relevel(dataset_tabs$negat_exp, ref = "1")
  
  
  # Fit satuarted model (XYZ)
  sat.model <- glm(Freq ~ varname * id_sex * negat_exp, 
                   family = poisson(), data = dataset_tabs, subset = Freq > 0)
  
  # Fit complete independence (X, Y, Z)
  
  compl.indep <- glm(Freq ~ varname + id_sex + negat_exp, 
                     family = poisson(), data = dataset_tabs, subset = Freq > 0)
  
  # Joint independence (X, YZ)
  joint.indep <- glm(Freq ~ varname + id_sex + negat_exp + 
                       negat_exp*varname, 
                     family = poisson(), data = dataset_tabs, subset = Freq > 0)
  
  # Conditional independence (XZ, YZ)
  condi.indep <- glm(Freq ~ varname + id_sex + negat_exp + 
                       negat_exp*varname + negat_exp*id_sex, 
                     family = poisson(), data = dataset_tabs, subset = Freq > 0)
  
  # Homogeneous Association (XZ, YZ, XY)
  
  homog.associ <- glm(Freq ~ (varname + id_sex + negat_exp)^2, 
                      family = poisson(), data = dataset_tabs, subset = Freq > 0)
  
  # Model list
  
  models_list <- list(sat.model, compl.indep, joint.indep, condi.indep, homog.associ)
  
  
  # Model fit
  
  models_fit <- vcdExtra::LRstats(glmlist(sat.model,
                                          compl.indep,
                                          joint.indep, 
                                          condi.indep,
                                          homog.associ))
  models_fit_dt <- as.data.frame(models_fit)
  
  if(!file.exists(here("outcomes", paste0("best_risk_model_fit_", i, ".docx")))) {
    
    models_fit_fx <- flextable(models_fit_dt)
    
    # create a word document and add the table
    doc_models_fit_fx <- read_docx() |> 
      body_add_flextable(value = models_fit_fx) 
    
    # Save the table 
    print(doc_models_fit_fx, target = here("outcomes", paste0("best_risk_model_fit_", i, ".docx")))
    
    
  }
  
  # Select the model based on AIC
  
  lowest_AIC_index <- which.min(models_fit$AIC)
  
  best_model <- models_list[[lowest_AIC_index]]
  
  best_model_typ_risk[[paste0("best_model_risk", i)]] <- best_model
  
}


transfor_result_risk <- list()

for (i in seq_along(best_model_typ_risk)) {
  
  transfor_result_risk[[i]] <- result_transfor_tb(i, contingency_tabs_risk, best_model_typ_risk)
}

transfor_result_risk_all <- bind_rows(transfor_result_risk)

# Check the unique values in type_of_victimization
actual_levels <- unique(transfor_result_risk_all$type_of_victimization)

# Changing order 
transfor_result_risk_all$type_of_victimization = factor(transfor_result_risk_all$type_of_victimization, 
                                                         levels = actual_levels)

# Pattern to remove
pattern_remove <- "cat_percv_risk_" 

# Rename interactions 
transfor_result_risk_all <- transfor_result_risk_all |> 
  mutate(worry_offense = gsub(pattern_remove, "", type_of_victimization),
         interactions = case_when(
           worry_offense == "cyber" &
             models_interaction == "varnameYes:id_sexMale" ~ "Perc. Risk Cyber: Perc. Risk x Male",
           worry_offense == "cyber" &
             models_interaction == "varnameYes:negat_exp0" ~ "Perc. Risk Cyber: Perc. Risk x No Negative Exp.",
           worry_offense == "cyber" &
             models_interaction == "id_sexMale:negat_exp0" ~ "Perc. Risk Cyber: Male x No Negative Exp.",
           worry_offense == "bully" &
             models_interaction == "varnameYes:id_sexMale" ~ "Perc. Risk Bullying: Perc. Risk x Male",
           worry_offense == "bully" &
             models_interaction == "varnameYes:negat_exp0" ~ "Perc. Risk Bullying: Perc. Risk x No Negative Exp.",
           worry_offense == "bully" &
             models_interaction == "id_sexMale:negat_exp0" ~ "Perc. Risk Bullying: Male x No Negative Exp.",
           worry_offense =="partner" &
             models_interaction == "varnameYes:id_sexMale:negat_exp0" ~ "Perc. Risk Dating Violence: Perc. Risk x Male x No Negative Exp.",
           worry_offense == "partner" &
             models_interaction == "varnameYes:id_sexMale" ~ "Perc. Risk Dating Violence: Perc. Risk x Male",
           worry_offense == "partner" &
             models_interaction == "varnameYes:negat_exp0" ~ "Perc. Risk Dating Violence: Perc. Risk x No Negative Exp.",
           worry_offense == "partner" &
             models_interaction == "id_sexMale:negat_exp0" ~ "Perc. Risk Dating Violence: Male x No Negative Exp.",
           worry_offense =="sex_vio" &
             models_interaction == "varnameYes:id_sexMale:negat_exp0" ~ "Perc. Risk Sexual Harassment: Perc. Risk x Male x No Negative Exp.",
           worry_offense == "sex_vio" &
             models_interaction == "varnameYes:id_sexMale" ~ "Perc. Risk Sexual Harassment: Perc. Risk x Male",
           worry_offense == "sex_vio" &
             models_interaction == "varnameYes:negat_exp0" ~ "Perc. Risk Sexual Harassment: Perc. Risk x No Negative Exp.",
           worry_offense == "sex_vio" &
             models_interaction == "id_sexMale:negat_exp0" ~ "Perc. Risk Sexual Harassment: Male x No Negative Exp."))

# Save the table as supplementary Table S3
if(!file.exists(here("outcomes", paste0("log_lin_res_risk.docx")))) {
  
  transfor_result_risk_fx <- flextable(transfor_result_risk_all)
  
  # create a word document and add the table
  doc_transfor_result_risk <- read_docx() |> 
    body_add_flextable(value = transfor_result_risk_fx) 
  
  # Save the table 
  print(doc_transfor_result_risk, target = here("outcomes", paste0("log_lin_res_risk.docx")))
  
  
}

transfor_result_risk_all <- transfor_result_risk_all |> 
  mutate(interactions_gen = extract_after_colon(interactions),
         perc_risk = extract_before_colon(interactions))


# Forestplot worry about points = interactions y =  models. 
forest_plot2 <- transfor_result_risk_all |>  
  ggplot(aes(y = perc_risk, x = cohens_d, color = interactions_gen, fill = interactions_gen)) +
  geom_point(shape = 21, size = 4, position = position_dodge(width = 0.75)) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high),
  position = position_dodge(width = 0.75), linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(y = "", 
       x = "Cohen's d",
       fill = "",
       color = "") +
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        text = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 2), 
         color = guide_legend(nrow = 2)) +
  scale_fill_viridis(discrete = TRUE,
                     option = "D",
                     begin = 0,
                     end = 0.9) +
  scale_color_viridis(discrete = TRUE,
                      option = "D",
                      begin = 0,
                      end = 0.9) 



# Figure 3
plot_grid(forest_plot1, forest_plot2, labels = "AUTO")

ggsave(here("submission", "TrinidadFigure3.tiff"), dpi = 300, 
       scale = 2, width = 10, height = 3)


sessionInfo()
