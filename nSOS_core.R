#' Next-Generation Specific Oxidative Stress Index (nSOS)
#' Core Computational Framework
#' License: MIT
#' Author: Independent Researcher, Sofia, Bulgaria

library(dplyr)
library(tidyr)

#' Calculate Non-Linear Fold-Change for a Single Biomarker
calculate_fold_change <- function(biomarker_vector, control_mean) {
  if (is.na(control_mean) || control_mean == 0) {
    stop("Critical Error: Control mean is invalid (NA or Zero). Check reference site data.")
  }
  
  ratio <- biomarker_vector / control_mean
  fc <- ifelse(ratio >= 1, ratio, -1 / ratio)
  return(abs(fc))
}

#' Compute the nSOS Index for a Multi-Biomarker Matrix
calculate_nsos <- function(data, biomarker_cols, control_site_name) {
  
  if (!"Site" %in% colnames(data)) {
    stop("Input data must contain a column named 'Site'.")
  }
  if (!control_site_name %in% data$Site) {
    stop(paste0("Specified control site '", control_site_name, "' not found in data."))
  }
  
  data <- data %>% mutate(across(all_of(biomarker_cols), ~as.numeric(as.character(.))))
  
  # Calculate reference baselines
  control_means <- data %>%
    filter(Site == control_site_name) %>%
    summarise(across(all_of(biomarker_cols), mean, na.rm = TRUE))
  
  fc_data <- data
  for (col in biomarker_cols) {
    fc_col_name <- paste0(col, "_FC")
    c_mean <- as.numeric(control_means[[col]])
    fc_data[[fc_col_name]] <- sapply(data[[col]], calculate_fold_change, control_mean = c_mean)
  }
  
  fc_cols <- paste0(biomarker_cols, "_FC")
  
  final_results <- fc_data %>%
    rowwise() %>%
    mutate(
      nSOS_Score = sqrt(mean(c_across(all_of(fc_cols))^2, na.rm = TRUE))
    ) %>%
    ungroup()
  
  return(final_results)
}

#' Generate Biologically Validated Mock Biomarker Data
generate_mock_data <- function() {
  set.seed(42)
  
  ref <- data.frame(
    Site = rep("Reference", 10),
    SOD   = rnorm(10, mean = 10, sd =  1),
    CAT   = rnorm(10, mean = 50, sd =  5),
    GST   = rnorm(10, mean = 5,  sd =  0.5),
    LPO   = rnorm(10, mean = 1.2, sd = 0.1)
  )
  
  site_a <- data.frame(
    Site = rep("Site_A", 10),
    SOD   = rnorm(10, mean = 25,  sd = 3),   
    CAT   = rnorm(10, mean = 120, sd = 12), 
    GST   = rnorm(10, mean = 14,  sd = 1.5), 
    LPO   = rnorm(10, mean = 3.5, sd = 0.4)  
  )
  
  site_b <- data.frame(
    Site = rep("Site_B", 10),
    SOD   = rnorm(10, mean = 2.5, sd = 0.3), 
    CAT   = rnorm(10, mean = 48,  sd = 4),   
    GST   = rnorm(10, mean = 20,  sd = 2),   
    LPO   = rnorm(10, mean = 4.8, sd = 0.5)  
  )
  
  mock_data <- rbind(ref, site_a, site_b)
  mock_data$ID <- 1:nrow(mock_data)
  
return(mock_data %>% select(ID, Site, SOD, CAT, GST, LPO))
  
  }