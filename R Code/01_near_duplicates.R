################################################################################
# Title: Examining Near Duplicates 
# Syntax source is based on the code by Jake Day, Jon Brauer, Mjaa Kotlja (2023, October 10)
# https://reluctantcriminologists.com/blog-posts/%5B8%5D/dup-series/duplicates
################################################################################
# Load packages
library(tidyverse)
library(here)
library(labelled)
library(janitor)
library(report)
library(extRemes)
library(qqplotr)

# Load new data 
mydata <- readRDS(here("data_ori", "name_data.rds"))

# Remove label

nolabel_data <- remove_val_labels(mydata) |> 
  mutate(Edad_Recodificada = case_when(
    Edad_Recodificada == 15 ~ 2,
    TRUE ~ Edad_Recodificada))

# Select only the data for worry about risk perceptions
data <- nolabel_data |> 
  select(2:7, 14, 40:57, 59:80, 83:103, 106:123, 126, 128:234)


# Checking exact duplication ----------------------------------------------

data |> 
  get_dupes(-c(CODIGO, Id.sexual))


# Checking near duplication -----------------------------------------------

# Data pre-processing
data_clean <- data |> 
  select(-c(CODIGO, Id.sexual)) |> 
  mutate(
    across(.cols = 1:192, ~na_if(., 999)))

# Identifying variables with more than 10% of missing as Kuriakose and Robbins (2016) 

missing_percentages <- map(data_clean, ~mean(is.na(.))) |> 
  as_tibble() |> 
  pivot_longer(
    cols = everything(),
    names_to = "items",
    values_to = "perc_missing") |> 
  filter(perc_missing >= .10) |> 
  arrange(desc(perc_missing))

# Select the names of the variables
missing_names <- missing_percentages[[1]]

# Removing those variables with more than 10 % of missing data
data_clean_missing_10 <- data_clean |> 
  select(-all_of(missing_names))

# Identify cases (rows) with more than 25 % missing  (K&R, 2016 criteria)
data_missing_cases <- map_dbl(seq_len(nrow(data_clean_missing_10)), 
                              ~ mean(is.na(data_clean_missing_10[., ]))) |> 
  as_tibble() |> 
  mutate(id = row_number()) |> 
  rename(perc_cascemiss = value) |> 
  filter(perc_cascemiss >= .25) |> 
  arrange(desc(perc_cascemiss))

# Extract the row_number with more than 25 % missing
missing_rows_num <- data_missing_cases[[2]]

# Remove rows with more than 25 % missing
data_clean_missing_row <- data_clean_missing_10 |> 
  mutate(id.rownum = row_number()) |> 
  filter(!id.rownum %in% missing_rows_num) |> 
  select(id.rownum, everything()) # Id first for the percent_match function 

# Maximum percent match
percent_match <- function (x) {
  #  x <- bl_t
  t <- as.tibble(t(x))
  t[t == ""] <- NA
  id <- x[[1]]
  nr <- nrow(t)
  nc <- ncol(t)
  m <- vector(mode = "numeric", length = nc-1)
  pm <- vector(mode = "numeric", length = nc)
  id_m <- vector(mode = "character", length = nc)
  
  for (i in 1:nc) {
    # i = 1
    # i = 4
    n <- t[,i]
    tt <- t[,-i]
    for (j in 1:(nc-1)) {
      m[j] <- sum(n[,1] == tt[,j], na.rm = TRUE) / sum(!is.na(n[,1]))
    }
    pm[i] <- max(m)
    id_m[i] <- tt[1,which.max(m)]
    if (i %% 10 == 0) { 
      print(i) 
    }
  }
  out <- tibble(id, pm, id_m)
  return(out)
}

# Percent matc
match_percentage <- percent_match(data_clean_missing_row)

# Identify those cases with a 85 % as in K&R (2016)
match_percentage_85 <- match_percentage |> 
  filter(pm >= .85) |> 
  arrange(desc(pm))

match_percentage_80 <- match_percentage |> 
  filter(pm >= .80) |> 
  arrange(desc(pm))

# Descriptive stats
descipt_pct_match <- match_percentage |> 
  select(pm) |> 
  report_table()  # Results present lower values than K & R simulaiton (mean = .66)


# Plotting the distribution of maximum proporiton

ggplot( data = match_percentage, aes(x = pm)) +
  geom_histogram(aes(y = after_stat(density)),
                 color = "white", fill = "darkred") +
  scale_x_continuous(limits = c(0.2, 0.99), breaks = seq(0, 1, .1)) +
  labs(title = "Distribution of Maximum Proportion Match",
       y = "Density",
       x = "Maximum Proportion Match") +
  theme_minimal() +
  theme(panel.grid = element_blank())
  
# Fitting Gumbel's distribution


fit_gumbel <- fevd(x = match_percentage$pm, type = "Gumbel",
     time.units = NULL, period.basis = NULL)  

# Plot Quantile Plot of Gumbel distribution. 

gumbel_parameter <- as_tibble_row(distill(fit_gumbel))

gumbel_qqplot <- ggplot(data = match_percentage, aes(sample = pm)) + 
  stat_qq_band(bandType = "boot", fill = "#218F8B", alpha = 0.25, identity = TRUE, B = 1000, 
               distribution = "evd",
               dparams = list(loc = gumbel_parameter$location, scale = gumbel_parameter$scale)) + 
  stat_qq_line(color = "#218F8B", identity = TRUE,
               distribution = "evd", 
               dparams = list(loc = gumbel_parameter$location, scale = gumbel_parameter$scale)) + 
  stat_qq_point(color = "#440154" , alpha = 0.25, 
                distribution = "evd", 
                dparams = list(loc = gumbel_parameter$location, scale = gumbel_parameter$scale)) +
  labs(title = "Quantile-Quantile Plot of Gumbel Distribution",
       x = "Gumbel Distribution Quantiles",
       y = "Maximum Proportion Match Quantiles") + 
  scale_y_continuous(breaks = seq(.2,1,.1)) +
  scale_x_continuous(breaks = seq(.2, 1, .2)) + 
  coord_cartesian(ylim = c(.2, 1), xlim = c(.2, 1)) + 
  theme_light() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())

gumbel_qqplot 

# Fitting an alternative extreme values distribution such as Generalized Extreme Values (GEV)

fit_gev <- fevd(x = match_percentage$pm, type = "GEV",
                      time.units = NULL, period.basis = NULL)  

gev_parameters <- as_tibble_row(distill(fit_gev))


gev_qqplot <- ggplot(data = match_percentage, aes(sample = pm)) + 
  stat_qq_band(bandType = "boot", fill = "#218F8B", alpha = 0.25, identity = TRUE, B = 1000, 
               distribution = "evd",
               dparams = list(loc = gev_parameters$location, scale = gev_parameters$scale, shape = gev_parameters$shape)) + 
  stat_qq_line(color = "#218F8B", identity = TRUE,
               distribution = "evd", 
               dparams = list(loc = gev_parameters$location, scale = gev_parameters$scale, shape = gev_parameters$shape)) + 
  stat_qq_point(color = "#440154" , alpha = 0.25, 
                distribution = "evd", 
                dparams = list(loc = gev_parameters$location, scale = gev_parameters$scale, shape = gev_parameters$shape)) + 
  labs(title = "Quantile-Quantile Plot of GEV Distribution",
       x = "Generalized Extreme Value Distribution Quantiles",
       y = "Maximum Proportion Match Quantiles") + 
  scale_y_continuous(breaks = seq(.2,1,.1)) +
  coord_cartesian(ylim = c(.2, 1) , xlim = c(0.2, 1)) + 
  theme_light() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())

gev_qqplot

# Simulating data an comparing

# Empirical Distribution Data
pctmatch <- match_percentage |> 
  select(pm)  |> 
  rename(pctmatch = pm) |> 
  mutate(dist = "Empirical")

# Gumbel Distribution Data
set.seed(7734)
gumbel <- revd(10000,loc = gumbel_parameter$location, scale = gumbel_parameter$scale) |> 
  as_tibble() |> 
  rename(pctmatch = value) |> 
  mutate(dist = "Gumbel")

# Weibull Distribution 
set.seed(7734)
weibull <- revd(10000,loc = gev_parameters$location, scale = gev_parameters$scale, shape = gev_parameters$shape) %>%
  as_tibble() |> 
  rename(pctmatch = value) |> 
  mutate(dist = "Weibull")

distcomp <- bind_rows(pctmatch, gumbel, weibull)

#Compare theoretical distributions to empirical
distcomp_plot <- ggplot(data = distcomp, aes(color = dist)) +
  geom_density(aes(x = pctmatch), linewidth = 1) + 
  # geom_vline(aes(xintercept = 0.70), color = "grey", linetype = 2, linewidth = .5) +
  # geom_vline(aes(xintercept = 0.80), color = "grey", linetype = 2, linewidth = .5) +
  scale_color_manual(values = c("#440154", "#FDE725FF", "#218F8B")) +
  labs(title = "Density Plot of Emperical Distribution and Theoretical Extreme Value Distributions",
       x = "Maximum Proportion Match",
       y = "Density",
       color = "Distribution") + 
  scale_x_continuous(breaks = seq(.2, 1.8, .2)) + 
  theme_light() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank())

distcomp_plot 
