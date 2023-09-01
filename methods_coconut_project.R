##############################################################################
# TESTING THE YIELD SENSITIVITY OF COCONUT PLANTS TO DROUGHT
##############################################################################

##############################################################################
# PROJECT GOALS
##############################################################################

# 1) DEFINE DROUGHT
# 2) DEFINE YIELD SENSITIVITY 
# 3) MEASURE SENSITIVITY TO DROUGHT
# 4) IDENTIFY SENSITIVE AND INSENSITIVE TREES
# 4.1) IDENTITY "SUPER PRODUCING" TREES

##############################################################################
rm(list=ls())
getwd()

setwd("Users/tomversluys/Documents/Coconut_data/Yield_data")

# load packages
##############################################################################
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(purrr)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(jtools)
library(tidyverse)
library(lubridate)
library(nlme)
library(lme4)
library(mgcv)
library(interactions)
library(SPEI)
library(tidymodels)
library(broom.mixed)
library(sjPlot)
library(glmmTMB)

## define functions
##############################################################################
rmse_mae_smape <- function(x){
  residuals <- residuals(x)
  mae <- mean(abs(residuals))
  rmse_value <- sqrt(mean(residuals^2))
  # Extracting predicted values
  predicted_values <- fitted(x)
  # Derive observed (true) values from predicted values and residuals
  observed_values <- predicted_values + residuals
  # Calculate sMAPE
  smape <- (2 * sum(abs(predicted_values - observed_values))) / 
    sum(predicted_values + observed_values) * 100
  return(list(RMSE = rmse_value, sMAPE = smape, MAE = mae))
}

rmse_mae_smape_logged <- function(x){
  # Get residuals and predicted values on the log scale
  residuals_log <- residuals(x)
  predicted_values_log <- fitted(x)
  # Back-transform to original scale
  predicted_values <- exp(predicted_values_log)
  observed_values <- exp(predicted_values_log + residuals_log)
  # Calculate MAE
  mae <- mean(abs(observed_values - predicted_values))
  # Calculate RMSE
  rmse_value <- sqrt(mean((observed_values - predicted_values)^2))
  # Add a small constant to the denominator to avoid potential divisions by zero
  epsilon <- 1e-10
  # Calculate sMAPE
  smape <- (2 * sum(abs(predicted_values - observed_values))) / 
    (sum(predicted_values + observed_values) + epsilon) * 100
  return(list(c(RMSE = rmse_value, sMAPE = smape, MAE = mae)))
}

get_avg_slope <- function(a_fixed, b_fixed, a_random, b_random, seq) {
  # Compute the combined coefficients
  a_total <- a_fixed + a_random
  b_total <- b_fixed + b_random
  # Compute the slope for each x in seq
  slopes <- 2 * a_total * seq + b_total
  # Return the average slope
  return(mean(slopes))
}

## read in coconut yield data
##############################################################################
df1_1986 <- fread("B022_1986_filtered.csv") %>% mutate(year_planted = 1986)
df2_1988 <- fread("B023_1988_filtered.csv") %>% mutate(year_planted = 1988)
df3_1992 <- fread("B084_1992_drop_month1_filtered.csv") %>% mutate(year_planted = 1992)
df4_1995 <- fread("B023_1995_drop_month11_filtered.csv") %>% mutate(year_planted = 1995)
df5_1997 <- fread("B052_1997_filtered.csv") %>% mutate(year_planted = 1997)
df6_1998 <- fread("B034_1998_filtered.csv") %>% mutate(year_planted = 1998)
df7_2002 <- fread("B050_2002_drop_month_9_11_filtered.csv") %>% mutate(year_planted = 2002)

# define function to select desired columns
##############################################################################
col_select_function <- function(df){
  df <- df %>% dplyr::select(tree_id, index, field, row, tree, year_planted, year, month,
                      nb_fr) %>% data.table()
}
colnames(df1_1986)
# subset data 
##############################################################################
df1_1986 <- col_select_function(df1_1986)
df2_1988 <- col_select_function(df2_1988)
df3_1992 <- col_select_function(df3_1992)
df4_1995 <- col_select_function(df4_1995)
df5_1997 <- col_select_function(df5_1997) 
df6_1998 <- col_select_function(df6_1998)
df7_2002 <- col_select_function(df7_2002)

# combine data
##############################################################################
df_join <- rbind(df1_1986, df2_1988, df3_1992, df4_1995, df5_1997, df6_1998, df7_2002)

# replace missing values with NAs
df_join[] <- lapply(df_join, function(x) ifelse(x == "", NA, x))

# examine data
##############################################################################
# note: original tree_id variable includes year of observation, not of planting - is this correct?
# create new tree ids to check accuracy of original "tree_id" variable
df_join1 <- df_join %>% mutate(new_tree_id = paste(year_planted, field, row, tree, sep = "_"))

# check number of unique trees defined by old variable
##############################################################################
length(unique(df_join1$tree_id)) # 41953 - is this more than expected?

# repeat with new variable 
##############################################################################
length(unique(df_join1$new_tree_id)) # 6585 - 86% disappear


## read in exclusions data
##############################################################################
exclusions <- fread("Sequenced_palms.csv") %>% mutate(exclusions_tree_id = paste(Year.planted, Row, Plant, sep = "_"))

# create matching ID for main data
df_join1 <- df_join1 %>% mutate(exclusions_tree_id = paste(year_planted, row, tree, sep = "_"))

length(unique(df_join1$exclusions_tree_id)) # 6581; this closely matches the number above

df_join1 <- df_join1 %>% filter(! exclusions_tree_id %in% exclusions$exclusions_tree_id) %>% 
  data.table()

length(unique(df_join1$exclusions_tree_id)) # 377 lost - missing just one

# explore the spatial structure of the data
##############################################################################
df_join1 %>% group_by(field) %>% summarise(row_no  = length(unique(row)),
                                               tree_no  = length(unique(new_tree_id)),
                                               year_planted_no  = length(unique(year_planted)), # year planted mainly indexes field
                                               year_planted  = unique(year_planted)) # but field 23 has trees planted in 1995 and 1998

df_join_plot %>% group_by(field, row) %>% summarise( # each field has rows, each with 26 trees
                                              tree_no  = length(unique(new_tree_id)),
                                              year_planted_no  = length(unique(year_planted)))

# check number of observations per tree
##############################################################################
df_join1 <- df_join1 %>% arrange(new_tree_id) %>% 
  group_by(new_tree_id, year, month) %>% # check only 1 observation per tree per year and month
  mutate(num = n()) %>% # some have too many, so remove
  filter(num == 1) %>% 
  dplyr::select(-num) %>%
  group_by(new_tree_id) %>% 
  mutate(observations_per_tree = n()) %>%
  filter(observations_per_tree > 1) %>% # remove trees with only a single observation
  data.table() 

df_join1 %>% group_by(observations_per_tree) %>% 
  summarise(l = length(unique(new_tree_id))) # between  8 and 84 obs per tree

# check distribution of yield overall
##############################################################################
df_join2 <- df_join1 %>% filter(! nb_fr == "dmf") %>% mutate(num_fruit = as.numeric(nb_fr)) %>% drop_na(num_fruit) %>% 
  dplyr::select(new_tree_id, year, month, year_planted, num_fruit, observations_per_tree, field, plot, row, tree) %>%
  data.table()

# histogram for overall yield
##############################################################################
(hist_overall_yeild <- df_join2 %>% group_by(new_tree_id) %>% slice(1) %>% 
   ggplot(aes(x = num_fruit)) + geom_histogram()) # roughly log normal

df_join2 <- df_join2 %>% mutate(m = mean(num_fruit, na.rm = T), # remove data points > 5sd away from the mean
                                sd = sd(num_fruit, na.rm = T)) %>%
  filter(! num_fruit > (m+(5*sd))) %>% dplyr::select(-c(m, sd))

# boxplot
##############################################################################
(boxplot_overall_yeild <- ggplot(df_join2, aes(x = num_fruit)) + 
   geom_boxplot()) # many outliers

df_join2 %>% group_by(field) %>%
  summarise(m = round(mean(num_fruit),0), v = round(sd(num_fruit),0))

# now look at variation within trees
##############################################################################
df_join2 <- df_join2 %>% group_by(new_tree_id) %>%
  mutate(m = round(mean(num_fruit),0), v = round(sd(num_fruit),0)) %>% data.table()

df_join2 <- df_join2 %>% arrange(new_tree_id, year)

# histogram for between tree variability
##############################################################################
(hist_between_tree_variability <- df_join2 %>% group_by(new_tree_id) %>% slice(1) %>%
   ggplot(aes(x = m)) + 
   geom_histogram()) # there's high variation in average yield

# histogram for within tree variability
##############################################################################
(hist_within_tree_variability <- df_join2 %>% group_by(new_tree_id) %>% slice(1) %>%
   ggplot(aes(x = v)) + 
   geom_histogram(bins = 20)) # some trees are far more variable in their yield than other

# variation by planting year
##############################################################################
df_join2 %>% mutate(tree_age = year - year_planted) %>% 
  group_by(tree_age) %>% 
  mutate(f = scale(num_fruit)) %>% 
  group_by(year_planted) %>%
  summarise(m = mean(f, na.rm = T)) # after adjusting for tree age, there's variation by planting year

# variation by observation year and month
##############################################################################
df_join2 %>% mutate(tree_age = year - year_planted) %>% 
  group_by(tree_age, year_planted) %>%
  mutate(f = scale(num_fruit)) %>% 
  # group_by(year) %>% # no variation by year of observation
  group_by(month) %>% # but variation by month of observation (late months have low yields)
  summarise(m = mean(f, na.rm = T)) 

# read in climate data
##############################################################################
# read in via loop as there are many files
getwd()
# setwd("Yield_data/Climate_data")
final_df <- data.frame()
# Assume the files are for the years 1990 to 2020
for(year in 1992:2021) {
  # set file names 
  filename <- paste0("Climatic_", year, ".csv")
  # read files
  data <- fread(filename)
  data <- data %>% dplyr::select(Month, Day, Pluie, Temperature, "Insolat.") %>% 
    data.table()

  names(data) <- c("month", "day", "rainfall", "temperature", "sunlight")

  data[] <- lapply(data, function(x) ifelse(x == "", NA, x))
  data <- data %>% 
    mutate(rainfall = as.numeric(rainfall), 
           rainfall = ifelse(is.na(rainfall), 0, rainfall)
           ) %>%
    drop_na(rainfall, sunlight, temperature) %>% 
    group_by(month) %>% 
    mutate(mean_monthly_rainfall = mean(rainfall, na.rm = T),
           total_monthly_rainfall = sum(rainfall, na.rm = T),
           mean_monthly_sunlight = mean(sunlight, na.rm = T),
           mean_monthly_temperature = mean(temperature, na.rm = T), 
           max_monthly_temperature = max(temperature, na.rm = T),
           min_monthly_temperature = min(temperature, na.rm = T)) %>%
    ungroup() %>%
    mutate(
      month = match(tolower(month), tolower(month.name)),
      year = year) %>%
    arrange(year, month) %>%
    data.frame()

  final_df <- rbind(final_df, data)
}

climate_data <- final_df %>% 
  group_by(year, month) %>%
  slice(1) %>%
  data.table()

# manipulate datasets to create climatic variables
##############################################################################
climate_annual <- climate_data %>% group_by(year) %>% 
  mutate(annual_temp = mean(mean_monthly_temperature, na.rm = T),
         total_annual_rainfall = sum(total_monthly_rainfall, na.rm = T)) %>%
  data.table()
lag <- climate_annual %>% group_by(year) %>% slice(1) %>% ungroup() %>%
  mutate(lag_total_annual_rainfall = lag(total_annual_rainfall, 1),
         lag2_total_annual_rainfall = lag(total_annual_rainfall, 2),
         lag3_total_annual_rainfall = lag(total_annual_rainfall, 3),
         lag4_total_annual_rainfall = lag(total_annual_rainfall, 4)) %>%
  dplyr::select(year, lag_total_annual_rainfall,lag2_total_annual_rainfall,
         lag3_total_annual_rainfall,lag4_total_annual_rainfall) %>%
  data.table()
climate_annual <- left_join(climate_annual, lag, by = "year")

## get PETs
##############################################################################
# https://www.rdocumentation.org/packages/SPEI/versions/1.7/topics/Potential%20evapotranspiration
# get climatic water balance (total precipitation - PET)
# https://github.com/sbegueria/SPEI
# https://drought.unl.edu/archive/Documents/NDMC/Workshops/136/Pres/Brian%20Fuchs--PDSI%20and%20scPDSI.pdf

climate_annual1 <- climate_annual %>% mutate(PET_hargreaves = hargreaves(min_monthly_temperature, max_monthly_temperature, 
                                                                         # Pre = total_monthly_rainfall, 
                                                                         lat = 3.23, na.rm = TRUE),
                                             PET_thornthwaite = thornthwaite(mean_monthly_temperature, lat = 3.23, na.rm = TRUE), 
                                             climatic_water_balance = total_monthly_rainfall - PET_hargreaves)

# here's how this works: 
##############################################################################
# e.g., take a 3 month sequence (jan, feb, march) of rainfall values and get its total 
# do this for all possible 3 month sequences, giving you X "totals"
# get the mean and sd of these totals 
# then for each total, turn it into a z score by substracting the mean of totals and dividing by the sd of totals
# see https://www.researchgate.net/figure/SPEI-drought-index-categories_tbl1_283244485 for "standard" interpretation

# this only accounts for rainfall
spi_6 = spi(climate_annual1$total_monthly_rainfall, 6)
spi_12 = spi(climate_annual1$total_monthly_rainfall, 12)
spi_18 = spi(climate_annual1$total_monthly_rainfall, 18)
spi_24 = spi(climate_annual1$total_monthly_rainfall, 24)

# this accounts for temperature
spei_6 = spei(climate_annual1$climatic_water_balance, 6)
spei_12 = spei(climate_annual1$climatic_water_balance, 12)
spei_18 = spei(climate_annual1$climatic_water_balance, 18)
spei_24 = spei(climate_annual1$climatic_water_balance, 24)

# add variables to df
climate_annual1 <- climate_annual1 %>% mutate(spi_6 = spi_6$fitted,
                                              spi_12 = spi_12$fitted,
                                              spi_18 = spi_18$fitted,
                                              spi_24 = spi_24$fitted,
                                              spei_6 = spei_6$fitted,
                                              spei_12 = spei_12$fitted,
                                              spei_18 = spei_18$fitted,
                                              spei_24 = spei_24$fitted)

# add two more custom variables (total water and signal to noise ratio, snr)
# snr accounts for total rainfall but also variance in rainfall
##############################################################################
compute_rolling_totals <- function(data, n) {
  totals <- numeric(nrow(data))
  for (i in n:nrow(data)) {
    totals[i] <- sum(data$total_monthly_rainfall[(i-n+1):i])
  }
  totals <- ifelse(totals == 0, NA, totals)
  total <- totals
  return(total)
}

climate_annual1$total_rainfall_6 <- compute_rolling_totals(climate_annual1, 6)
climate_annual1$total_rainfall_12 <- compute_rolling_totals(climate_annual1, 12)
climate_annual1$total_rainfall_18 <- compute_rolling_totals(climate_annual1, 18)
climate_annual1$total_rainfall_24 <- compute_rolling_totals(climate_annual1, 24)

compute_rolling_totals <- function(data, n) {
  totals <- numeric(nrow(data))
  for (i in n:nrow(data)) {
    totals[i] <- mean(data$mean_monthly_temperature[(i-n+1):i])
  }
  totals <- ifelse(totals == 0, NA, totals)
  total <- totals
  return(total)
}

climate_annual1$mean_temp_6 <- compute_rolling_totals(climate_annual1, 6)
climate_annual1$mean_temp_12 <- compute_rolling_totals(climate_annual1, 12)
climate_annual1$mean_temp_18 <- compute_rolling_totals(climate_annual1, 18)
climate_annual1$mean_temp_24 <- compute_rolling_totals(climate_annual1, 24)

compute_rolling_signal_noise_ratio <- function(data, n) {
  snr <- numeric(nrow(data))
  for (i in n:nrow(data)) {
    snr[i] <- sum(data$total_monthly_rainfall[(i-n+1):i]) /
      sd(data$total_monthly_rainfall[(i-n+1):i])
  }
  snr <- ifelse(snr == 0, NA, snr)
  snr <- snr
  return(snr)
}

climate_annual1$snr_6 <- compute_rolling_signal_noise_ratio(climate_annual1, 6)
climate_annual1$snr_12 <- compute_rolling_signal_noise_ratio(climate_annual1, 12)
climate_annual1$snr_18 <- compute_rolling_signal_noise_ratio(climate_annual1, 18)
climate_annual1$snr_24 <- compute_rolling_signal_noise_ratio(climate_annual1, 24)

compute_rolling_counts <- function(data, n) {
  counts <- numeric(nrow(data))
  for (i in n:nrow(data)) {
    window <- data$total_monthly_rainfall[(i-n+1):i]
    counts[i] <- sum(window < 150)
  }
  
  return(counts)
}

climate_annual1$months_under_100_6 <- compute_rolling_counts(climate_annual1, 6)
climate_annual1$months_under_100_12 <- compute_rolling_counts(climate_annual1, 12)
climate_annual1$months_under_100_18 <- compute_rolling_counts(climate_annual1, 18)
climate_annual1$months_under_100_24 <- compute_rolling_counts(climate_annual1, 24)

## compute rdi
climate_annual1$pet_6 <- zoo::rollapply(climate_annual1$PET_thornthwaite,
                                                     width=6, FUN=sum, align="right", fill=NA)
climate_annual1$pet_12 <- zoo::rollapply(climate_annual1$PET_thornthwaite,
                                        width=12, FUN=sum, align="right", fill=NA)
climate_annual1$pet_18 <- zoo::rollapply(climate_annual1$PET_thornthwaite,
                                         width=18, FUN=sum, align="right", fill=NA)
climate_annual1$pet_24 <- zoo::rollapply(climate_annual1$PET_thornthwaite,
                                         width=24, FUN=sum, align="right", fill=NA)

climate_annual1 <- climate_annual1 %>% mutate(ratio_6 = total_rainfall_6/pet_6, 
                                              rdi_6 = (ratio_6 - mean(ratio_6, na.rm = TRUE))/sd(ratio_6, na.rm = TRUE), 
                                              ratio_12 = total_rainfall_12/pet_12, 
                                              rdi_12 = (ratio_12 - mean(ratio_12, na.rm = TRUE))/sd(ratio_12, na.rm = TRUE),
                                              ratio_18 = total_rainfall_18/pet_18, 
                                              rdi_18 = (ratio_18 - mean(ratio_18, na.rm = TRUE))/sd(ratio_18, na.rm = TRUE),
                                              ratio_24 = total_rainfall_24/pet_24, 
                                              rdi_24 = (ratio_24 - mean(ratio_24, na.rm = TRUE))/sd(ratio_18, na.rm = TRUE))

# check distribution of new variables
##############################################################################
climate_annual1 %>% pivot_longer(cols = c(spi_6, spi_12, spi_18, spi_24), 
                                         names_to = "vars", values_to = "vals") %>%
  ggplot(aes(vals)) + geom_histogram() + theme_few() + facet_grid(~vars)

climate_annual1 %>% pivot_longer(cols = c(spei_6, spei_12, spei_18, spei_24), 
                                         names_to = "vars", values_to = "vals") %>%
  ggplot(aes(vals)) + geom_histogram() + theme_few() + facet_grid(~vars)

climate_annual1 %>% pivot_longer(cols = c(total_rainfall_6, total_rainfall_12, 
                                          total_rainfall_18, total_rainfall_24), 
                                 names_to = "vars", values_to = "vals") %>%
  ggplot(aes(vals)) + geom_histogram() + theme_few() + facet_grid(~vars)

climate_annual1 %>% pivot_longer(cols = c(snr_6, snr_12, 
                                          snr_18, snr_24), 
                                 names_to = "vars", values_to = "vals") %>%
  ggplot(aes(vals)) + geom_histogram() + theme_few() + facet_grid(~vars)

climate_annual1 %>% pivot_longer(cols = c(rdi_6, rdi_12, 
                                          rdi_18), 
                                 names_to = "vars", values_to = "vals") %>%
  ggplot(aes(vals)) + geom_histogram() + theme_few() + facet_grid(~vars)

# check missingness across months
##############################################################################
climate_annual1 %>% pivot_longer(cols = c(spi_6, spi_12, spi_18, spi_24), 
                                 names_to = "vars", values_to = "vals") %>%
  mutate(na_yes = as.factor(ifelse(is.na(vals), "na", "present"))) %>%
  group_by(vars, month, na_yes) %>%
  summarise(number_na = n()) %>%
  ggplot(aes(x = na_yes, y = number_na, fill = na_yes)) + 
  geom_bar(stat = "identity") + 
  facet_grid(vars ~ month)

climate_annual1 %>% pivot_longer(cols = c(spei_6, spei_12, spei_18, spei_24), 
                                 names_to = "vars", values_to = "vals") %>%
  mutate(na_yes = as.factor(ifelse(is.na(vals), "na", "present"))) %>%
  group_by(vars, month, na_yes) %>%
  summarise(number_na = n()) %>%
  ggplot(aes(x = na_yes, y = number_na, fill = na_yes)) + 
  geom_bar(stat = "identity") + 
  facet_grid(vars ~ month)

# no NAs for other variables

# join climate and yield datasets
##############################################################################
all_data <- inner_join(df_join2, climate_annual1, by = c("year", "month")) 

# create model dataset
##############################################################################
all_data_annual_sample <- all_data %>% 
  mutate(new_tree_id = as.factor(new_tree_id),
    num_fruit_log = as.numeric(log(num_fruit + 1)),
    tree_age = as.numeric(year - year_planted), 
         year_planted = as.factor(year_planted), # change var format for mixed models
         year = as.factor(year),
        field = as.factor(field),
    plot = as.factor(plot),
         month = as.factor(month)) %>%
  group_by(new_tree_id) %>% 
  mutate(rainfall_threshold_18 = ifelse(total_rainfall_18 < 2250, 1, 0),
         rainfall_threshold_12 = ifelse(total_rainfall_12 < 1500, 1, 0),
         crude_rainfall_threshold = ifelse(total_annual_rainfall < 1500, 1, 0),
         lag_crude_rainfall_threshold = ifelse(lag_total_annual_rainfall < 1500, 1, 0),
         lag2_crude_rainfall_threshold = ifelse(lag2_total_annual_rainfall < 1500, 1, 0),
         lag3_crude_rainfall_threshold = ifelse(lag3_total_annual_rainfall < 1500, 1, 0),
         lag4_crude_rainfall_threshold = ifelse(lag4_total_annual_rainfall < 1500, 1, 0)) %>%
  mutate(n = n()) %>%
  mutate(time = 1:length(new_tree_id)) %>%
  ungroup() %>%
  dplyr::select(new_tree_id, 
         tree_age_z, 
         field,
         row, tree,
         field, tree_age, year_planted, num_fruit, num_fruit_log, n, month, time, year,
         spi_6, spi_12, spi_18,  rdi_6, rdi_12, rdi_18, 
         spei_6, spei_12, spei_18, mean_temp_6, mean_temp_12, mean_temp_18, 
         total_rainfall_6, total_rainfall_12, total_rainfall_18, 
         snr_6, snr_12, snr_18, 
         threshold_rdi_12, threshold_rdi_18, threshold_spi_12, threshold_spi_18,
         rainfall_threshold_18, rainfall_threshold_12, threshold_spei_12, threshold_spei_18) %>%
  data.table()

length(unique(all_data_annual_sample$new_tree_id))

# get annualised data
##############################################################################
annualised <- all_data_annual_sample %>% pivot_longer(cols = c(spi_6, spi_12, spi_18, 
                                                                # spi_24, rdi_24, snr_24, spei_24, total_rainfall_24,mean_temp_24,
                                                               rdi_6, rdi_12, rdi_18, 
                                                                 spei_6, spei_12, spei_18, 
                                                               total_rainfall_6, total_rainfall_12, total_rainfall_18, 
                                                               mean_temp_6, mean_temp_12, mean_temp_18, 
                                                               snr_6, snr_12, snr_18), 
                                                                 names_to = "vars", 
                                                                 values_to = "vals") %>%
  group_by(new_tree_id, year, vars) %>% mutate(vals = mean(vals, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = vars, values_from = vals) %>%
  group_by(new_tree_id, year) %>%
  mutate(yearly_nuts = mean(num_fruit, na.rm = T), 
         yearly_nuts_log = mean(num_fruit_log),
         yearly_nuts_log_1 = log(yearly_nuts + 1),
         yearly_nuts_count = as.integer(round(yearly_nuts, 0))) %>%
  slice(1) %>% 
  group_by(new_tree_id) %>%
  mutate(time = 1:length(new_tree_id), 
         n = n()) %>%
  filter(n >= 12) %>%
  drop_na() %>%
  data.table()

# check that calculations manually using a random tree
annualised1 %>% filter(new_tree_id == "1986_22_18_19", year == "1994") %>%
  summarise(mean(yearly_nuts)) 
all_data_annual_sample %>% filter(new_tree_id == "1986_22_18_19", year == "1994") %>%
  summarise(mean(num_fruit)) 

# check number of trees
length(unique(annualised$new_tree_id))

# check yield distribution
hist(annualised$yearly_nuts_count)

# identify "super producers" in the top 25% for mean and bottom 50% for variance
##############################################################################

# approach 1: partial out effects of tree age, field, and collection year
##############################################################################
system.time(super_mod <- glmmTMB(yearly_nuts_count ~ poly(tree_age, 2, raw = TRUE) + field + as.factor(year) + 
                                            (1 | new_tree_id),
                                            family="nbinom2",
                                            control=glmmTMBControl(parallel=8),
                                            data=annualised1))
summary(super_mod)
p_residuals <- residuals(super_mod, type="pearson")
super_tree_approach_1 <- annualised1 %>% select(new_tree_id, yearly_nuts_count) %>% data.table() %>%
  mutate(p_resids = p_residuals) %>%
  group_by(new_tree_id) %>%  mutate(yearly_mean_resid_p = mean(p_resids, na.rm = T), 
                                    yearly_sd_resid_p = sd(p_resids, na.rm = T)) %>%
  group_by(new_tree_id) %>% slice(1) %>%
  data.table()

quantile(super_tree_approach_1$yearly_mean_resid_p, 0.8)
quantile(super_tree_approach_1$yearly_sd_resid_p, 0.5)
hist(super_tree_approach_1$yearly_mean_resid_p)


# approach 1: get z score by tree age, field, and collection year
##############################################################################
super_tree_approach_2 <- annualised1 %>% group_by(field, tree_age, year) %>%
  mutate(nuts_z = as.numeric(scale(yearly_nuts))) %>% 
  group_by(new_tree_id) %>%
  mutate(yearly_mean_nuts = mean(nuts_z, na.rm = T), 
         yearly_sd_nuts = sd(nuts_z, na.rm = T)) %>%
  ungroup() %>%
  group_by(new_tree_id) %>% 
  slice(1) %>%
  data.table()

quantile(super_tree_approach_2$yearly_mean_nuts, 0.8)
quantile(super_tree_approach_2$yearly_sd_nuts, 0.5)
hist(super_tree_approach_2$yearly_mean_nuts)

# using second approach, filter to identify super trees
##############################################################################
super1 <- super_tree_approach_2 %>%  
  mutate(super_producer = ifelse(yearly_mean_nuts > 0.5773538 & yearly_sd_nuts < 0.6231901, "super_good", "normal")) %>%
  filter(super_producer == "super_good") %>%
  group_by(new_tree_id) %>%
  slice(1) %>%
  data.table()

# plot super trees vs. normal trees to visualise effects
##############################################################################
super_plot <- super1 %>% 
  filter(super_producer == "super_good") %>%
  filter(new_tree_id %in% sample(unique(new_tree_id), 10)) %>%
  ggplot(aes(x = year, y = nuts_z, colour = new_tree_id)) + 
  geom_point() + 
  geom_line(aes(group = new_tree_id)) +  
  ylim(-3, 4) +
  theme_minimal() +
  theme(legend.position = "none")

normal_plot <- super1 %>% 
  filter(! super_producer == "super_good") %>%
  filter(new_tree_id %in% sample(unique(new_tree_id), 10)) %>%
  ggplot(aes(x = year, y = nuts_z, colour = new_tree_id)) + 
  geom_point() + 
  geom_line(aes(group = new_tree_id)) +  
  ylim(-3, 4) +
  theme_minimal() +
  theme(legend.position = "none")

ggarrange(super_plot, normal_plot)

# next, explore relationship between drought and yield
##############################################################################

# first, fit linear models for all variables and compare model fits
##############################################################################

model_function_linear <- function(data) {

  model_list <- c("spi_6", "spi_12",
                  "spi_18", "spi_24",
                  "rdi_6", "rdi_12", "rdi_18", "rdi_24",
                  "snr_6", "snr_12", "snr_18", "snr_24",
                  "spei_6", "spei_12", "spei_18", "spei_24",
                  "total_rainfall_6", "total_rainfall_12", "total_rainfall_18", "total_rainfall_24")
  coef_linear <- numeric(length(model_list))
  rmse_vals <- numeric(length(model_list))
  mae_vals <- numeric(length(model_list))
  smape_vals <- numeric(length(model_list))
  
  for (i in 1:length(model_list)) {
    try({
      formula_str <- paste("yearly_nuts_log_1 ~ poly(", model_list[i], ", 2, raw = TRUE) + poly(tree_age, 2, raw = TRUE) + field")
      formula <- as.formula(formula_str)
      random_str <- paste("~ poly(", model_list[i], ", 2, raw = TRUE) | new_tree_id")
      random_formula <- as.formula(random_str)
      model <- lme(formula,
                   random = random_formula,
                   method = "REML",
                   control = lmeControl(opt = "optim"),
                   correlation = corCAR1(form = ~ time | new_tree_id),
                   data = data)
      
      coef_linear[i] <- fixef(model)[2]
      rmse_vals[i] <- as.numeric(rmse_mae_smape(model)[1])
      mae_vals[i] <- as.numeric(rmse_mae_smape(model)[2])
      smape_vals[i] <- as.numeric(rmse_mae_smape(model)[3])
      
    }, silent = TRUE)
  }
  results <- data.frame(variable = model_list, coef_linear = coef_linear,
                        rmse = rmse_vals, mae = mae_vals, smape = smape_vals)
  return(results)
}

system.time(results_linear_mods <- model_function_linear(annualised1))


# first, fit quadratric models for all variables and compare model fits
##############################################################################
model_function_quadratic <- function(data) {
  
  model_list <- c("spi_6", "spi_12", "spi_18", "spi_24",
                  "rdi_6", "rdi_12", "rdi_18", "rdi_24",
                  "snr_6", "snr_12", "snr_18", "snr_24",
                  "spei_6", "spei_12", "spei_18", "spei_24",
                  "total_rainfall_6", "total_rainfall_12", "total_rainfall_18", "total_rainfall_24")
  
  coef_linear <- rep(NA, length(model_list))
  coef_quadratic <- rep(NA, length(model_list))
  rmse_vals <- rep(NA, length(model_list))
  mae_vals <- rep(NA, length(model_list))
  smape_vals <- rep(NA, length(model_list))
  
  for (i in 1:length(model_list)) {
    try({
      formula_str <- paste0("yearly_nuts_count ~ poly(", model_list[i], ", 2, raw = TRUE) + poly(tree_age, 2, raw = TRUE) + field + (poly(", model_list[i], ", 2, raw = TRUE) | new_tree_id)")
      formula_obj <- as.formula(formula_str)
      model <- glmmTMB(formula_obj, 
                       family = "nbinom2",
                       control = glmmTMBControl(parallel = 8),
                       data = data)
      
      coef_linear[i] <- fixef(model)$cond[2]
      coef_quadratic[i] <- fixef(model)$cond[3]
      rmse_vals[i] <- as.numeric(rmse_mae_smape(model)[1])
      mae_vals[i] <- as.numeric(rmse_mae_smape(model)[2])
      smape_vals[i] <- as.numeric(rmse_mae_smape(model)[3])
      
    }, silent = TRUE)
  }
  results <- data.frame(variable = model_list, coef_linear = coef_linear, coef_quadratic = coef_quadratic, 
                        rmse = rmse_vals, mae = mae_vals, smape = smape_vals)
  return(results)
}

system.time(results_quadratic_mods <- model_function_quadratic(annualised1))

# now examine models in greater depth
##############################################################################

# can also use log model if desired (below) 

# system.time(rdi_18_log <- lme(yearly_nuts_log_1 ~ poly(rdi_18, 2, raw = TRUE) + poly(tree_age, 2, raw = TRUE) +
#                                 random = ~ poly(rdi_18, 2, raw = TRUE) | new_tree_id),
#                                 method = "REML",
#                                 control = lmeControl(opt = "optim"),
#                                 correlation = corCAR1(form = ~ time | new_tree_id),
#                                 data=annualised1)

# first, rdi_18
system.time(rdi_18_nb <- glmmTMB(yearly_nuts_count ~ poly(rdi_18, 2, raw = TRUE) + poly(tree_age, 2, raw = TRUE) + 
                                             (poly(rdi_18, 2, raw = TRUE) | new_tree_id) + (poly(rdi_18, 2, raw = TRUE) | field), 
                                           family="nbinom2",
                                           control=glmmTMBControl(parallel=8),
                                           data=annualised1))
summary(rdi_18_nb)
plot_model(rdi_18_nb, type = "pred", terms = "rdi_18[all]")
rmse_mae_smape(rdi_18_nb)

# second, snr_18
system.time(snr_18_nb <- glmmTMB(yearly_nuts_count ~ poly(snr_18, 2, raw = TRUE) + poly(tree_age, 2, raw = TRUE) + field +
                                        (poly(snr_18, 2, raw = TRUE) | new_tree_id),
                                      family="nbinom2",
                                      control=glmmTMBControl(parallel=8),
                                      data=annualised1))
summary(snr_18_nb)
plot_model(snr_18_nb, type = "pred", terms = "snr_18[all]")
rmse_mae_smape(snr_18_nb)

# check out of sample accuracy 
##############################################################################

# check if slopes are similar for two drought metrics
##############################################################################
s1 <- tidy(rdi_18_nb, effects = "ran_vals") %>% filter(group == "new_tree_id", term == "poly(rdi_18, 2, raw = TRUE)2") %>%
  mutate(rdi_slope = round(estimate, 3)) 
s2 <- tidy(snr_18_nb, effects = "ran_vals") %>% filter(group == "new_tree_id", term == "poly(snr_18, 2, raw = TRUE)2") %>%
  mutate(snr_slope = round(estimate, 3)) 

cbind(s1, s2) %>% cor()

# check if signs are similar
cbind(s1, s2) %>% mutate(match = ifelse(c(rdi_count_plot > 0 & rdi_count > 0) | 
                                        c(rdi_count_plot < 0 & rdi_count < 0), 
                                        "match", "not")) %>%
  group_by(match) %>% summarise(n = n()) %>% ungroup() %>% 
  summarise(n, perc = n/sum(n))

cbind(s1, s2) %>% pivot_longer(cols = c("rdi_slope",  "snr_slope"),
                         values_to = "val", names_to = "model") %>%
  ggplot(aes(x = val)) +
  geom_histogram() + 
  facet_grid(~ model, scales = "free")

# also compare to zero inflated model
##############################################################################
system.time(rdi_18_nb_zi <- glmmTMB(yearly_nuts_count ~ poly(rdi_18, 2, raw = TRUE) + poly(tree_age, 2, raw = TRUE) + field +
                                         (poly(rdi_18, 2, raw = TRUE) | new_tree_id),
                                       ziformula = ~ (poly(rdi_18, 2, raw = TRUE)) + poly(tree_age, 2, raw = TRUE) + field, 
                                       family="nbinom2",
                                       control=glmmTMBControl(parallel=8),
                                       data=annualised1))

# fit's worse on zero inflated
rmse_mae_smape(rdi_18_nb)
rmse_mae_smape(rdi_18_nb_zi)

# check model diagnostics 
##############################################################################
library(DHARMa)
res_linear <- simulateResiduals(rdi_18_nb_test)
par(mar = c(1, 1, 1, 1))
plot(res_linear)
testResiduals(res_linear)
testZeroInflation(res_linear)

# Compute standardized residuals
std_res <- residuals(rdi_poly_nb1, type = "pearson")

# Flag observations with large residuals, for instance, beyond +/- 2.5 or 3
outliers <- which(abs(std_res) > 2.5)

# computer sensitivity 
##############################################################################
get_avg_slope <- function(a_fixed, b_fixed, a_random, b_random, x_range) {
  a = a_fixed + a_random
  b = b_fixed + b_random
  
  # For each value in the range, compute the slope and then take the average
  avg_slope <- mean(sapply(x_range, function(x) 2*a*x + b))
  return(avg_slope)
}

# get model fixed and random effects
##############################################################################
fixed_effects <- fixef(rdi_18_nb)
ranef_tree <- as.data.frame(ranef(rdi_18_nb))
ranef_rdi <- ranef_tree %>% 
  mutate(level = grp, estimate = as.numeric(condval)) %>%
  filter(term %in% c("poly(rdi_18, 2, raw = TRUE)1", "poly(rdi_18, 2, raw = TRUE)2"))

# get tree sensitivities
##############################################################################
tree_sensitivities_df <- ranef_rdi %>%
  group_by(level) %>%
  dplyr::summarise(
    a_random = sum(estimate[term == "poly(rdi_18, 2, raw = TRUE)2"]),
    b_random = sum(estimate[term == "poly(rdi_18, 2, raw = TRUE)1"]),
    linear_term = (b_random + fixed_effects[["cond"]][["poly(rdi_18, 2, raw = TRUE)1"]]),
    quadratic_term = (a_random + fixed_effects[["cond"]][["poly(rdi_18, 2, raw = TRUE)2"]]),
    avg_slope_dry = get_avg_slope(
      a_fixed = fixed_effects[["cond"]][["poly(rdi_18, 2, raw = TRUE)2"]],
      b_fixed = fixed_effects[["cond"]][["poly(rdi_18, 2, raw = TRUE)1"]],
      a_random, 
      b_random,
      seq(-1, 0, by = 0.01)  # dry range
    ),
    avg_slope_wet = get_avg_slope(
      a_fixed = fixed_effects[["cond"]][["poly(rdi_18, 2, raw = TRUE)2"]],
      b_fixed = fixed_effects[["cond"]][["poly(rdi_18, 2, raw = TRUE)1"]],
      a_random, 
      b_random, 
      seq(0, 1.5, by = 0.01)   # wet range
    )
  )

tree_sensitivities_df <- tree_sensitivities_df %>% mutate(new_tree_id = level)
# tree_sensitivities_df <- tree_sensitivities_df %>% 
#   mutate(new_tree_id = tree_id)

set.seed(123)  # Ensures reproducibility
sample_data <- annualised1 %>% 
  data.table()

predicted_values <- predict(rdi_18_nb, newdata = sample_data, type = "response")
sample_data$predicted <- predicted_values


# Helper function to count sequences of zeros
count_zero_sequences <- function(values) {
  rle_values <- rle(values == 0)
  sum(rle_values$lengths[rle_values$values] > 5)
}

sample_data1 <- sample_data %>% 
  left_join(tree_sensitivities_df, by = "new_tree_id") %>%
  group_by(new_tree_id) %>%
  mutate(yearly_mean_nuts = mean(yearly_nuts, na.rm = T), 
         yearly_sd_nuts = sd(yearly_nuts, na.rm = T)) %>%
  mutate(num_zeros = sum(yearly_nuts_count == 0),
         n = n(),
         zero_perc = num_zeros/n,
         zero_sequences = count_zero_sequences(yearly_nuts_count)) %>%
  filter(zero_sequences == 0) %>% # lots of years entirely zero means a problem
  ungroup() %>%
  mutate(rank_dry = dense_rank(-avg_slope_dry),
         rank_wet = dense_rank(avg_slope_wet)) %>%
  group_by(new_tree_id) %>%
  # filter(rank_dry < 50) %>% # trees respond badly to little water
  # filter(rank_dry > 3750) %>% # trees respond well to little water
  # filter(rank_wet < 100) %>% # low (~100) trees the respond badly to excess water
  # filter(rank_wet > 3750) %>% # high (3750) trees that respond well to excess water
  # filter(between(avg_slope_dry, -0.01, 0.01)) %>% # trees that are unresponsive to drought
  # filter(between(avg_slope_wet, -0.02, 0.02)) %>% # trees that are unresponsive to excess water
  
  # filter(quadratic_term < 0) %>%
  # filter(between(quadratic_term, -0.05, 0.05)) %>%
  # filter(between(avg_slope_dry, -0.05, 0.05)) %>%
  # filter(between(avg_slope_wet, -0.05, 0.05)) %>%
# filter(between(linear_term, -0.01, 0.01)) %>%
# filter(between(quadratic_term, -0.01, 0.01)) %>%
# filter(between(avg_slope_dry, -0.2, -0.05)) %>%
# filter(between(avg_slope_wet, -0.1, 0.1)) %>%
# filter(avg_slope_wet < -0.2) %>%
# filter(avg_slope_dry > -0.2) %>%
ungroup() %>%
  # filter(new_tree_id %in% sample(unique(new_tree_id), 20)) %>%
  dplyr::select(new_tree_id, rdi_18, yearly_nuts_log, yearly_nuts_count, predicted, year,
                linear_term, quadratic_term, rank_dry, rank_wet, avg_slope_dry, avg_slope_wet,
                # year_planted, field, row, tree, 
                yearly_mean_nuts, yearly_sd_nuts) %>%
  data.table() 

# visualise sensitivities
##############################################################################
ggplot(sample_data1, aes(x=rdi_18)) +
  geom_point(aes(y=yearly_nuts_count, alpha = 0.2), color="red", size=2) +   # actual data points in red
  geom_smooth(aes(y=predicted), color="blue", method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  # geom_smooth(aes(y=yearly_nuts_count), color="red", method = "gam", formula = y ~ s(x, bs = "cs"),se = FALSE) +
  geom_point(aes(y=predicted), color="blue", size=1) +# predicted polynomial trend in blue
  facet_wrap(~ new_tree_id) +
  theme_few() +
  labs(title="Actual vs Predicted for Sampled Trees",
       x="RDI Value",
       y="Yearly Nut Count (Log)") +
  geom_vline(xintercept = 0)
