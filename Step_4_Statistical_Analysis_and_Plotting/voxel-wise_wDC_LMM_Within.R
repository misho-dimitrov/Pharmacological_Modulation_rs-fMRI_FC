# Load required libraries
library(tidyverse)
library(neurobase)
library(readr)
library(lme4)
#library(ptestr)
library(pTFCE)

###########################
setwd('~/Downloads/PhD/Analysis/Tianeptine/Static_FC/Jupyter_test_notebook/DC_data/Whole_brain_wDC_Z/')

subj_list_raw <- list.files('~/Downloads/PhD/Analysis/Tianeptine/Static_FC/Jupyter_test_notebook/DC_data/Whole_brain_wDC_Z/')

# Discard unnecessary files and sessions that did not pass the QC
subj_list_raw <- purrr::discard(subj_list_raw,.p = ~stringr::str_detect(.x,"final_resampled_gm.nii.gz"))
subj_list_raw <- purrr::discard(subj_list_raw,.p = ~stringr::str_detect(.x,"td_asd_intersected_constrained_0.8.nii.gz"))
subj_list_raw <- purrr::discard(subj_list_raw,.p = ~stringr::str_detect(.x,"asd_td_intersected.nii.gz"))
subj_list_raw <- purrr::discard(subj_list_raw,.p = ~stringr::str_detect(.x,"BRCTRADA005C_wDC_Z.nii"))
subj_list_raw <- purrr::discard(subj_list_raw,.p = ~stringr::str_detect(.x,"BRCTRADA006B_wDC_Z.nii"))
subj_list_raw <- purrr::discard(subj_list_raw,.p = ~stringr::str_detect(.x,"BRCTRADA006C_wDC_Z.nii"))
subj_list_raw <- purrr::discard(subj_list_raw,.p = ~stringr::str_detect(.x,"BRCTRADA108D_wDC_Z.nii"))
subj_list_raw <- purrr::discard(subj_list_raw,.p = ~stringr::str_detect(.x,"BRCTRADA111B_wDC_Z.nii"))
subj_list_raw <- purrr::discard(subj_list_raw,.p = ~stringr::str_detect(.x,"BRCTRADA115C_wDC_Z.nii"))

# Use grep to find indices of elements that match the pattern
# Either NT or A!!!!!!!!!
# ........

pattern <- "^BRCTRADA1"  # Regular expression pattern
matching_indices <- grep(pattern, subj_list_raw)

# Subset subj_list based on the matching indices
subj_list <- subj_list_raw[matching_indices]

# Print the new list
print(subj_list)


subj_id_list <- gsub('.{10}$', '', subj_list)


###########################
# Preallocate memory for empty_df
empty_df <- matrix(NA, nrow = 0, ncol = 196629)
col_names <- sprintf("voxel[%s]", seq_len(196629))
colnames(empty_df) <- col_names

# Initialize a list to store data frames
df_list <- list()

# Loop over subj_id_list
for (i in seq_along(subj_list)) {
  # Read the nifti file
  wDC <- readnii(subj_list[i])
  
  # Convert to data frame and transpose
  wDC_df_T <- as.data.frame(t(data.frame(wDC = c(wDC))))
  
  # Set row names
  row.names(wDC_df_T) <- subj_id_list[i]
  
  # Store data frame in the list
  df_list[[i]] <- wDC_df_T
  
  # Clear unnecessary variables
  rm(wDC, wDC_df_T)
}

# Combine all data frames in the list
final_df <- do.call(rbind, df_list)

# Clear unnecessary variables
rm(df_list)

# Now you have the final_df with the combined data BUT without predictors and covariates
# Add them
# Read the IV/covariate dataframe from a CSV file
# Either NT or A!!!!!!!!!
# ........
df_iv_cov <- read.csv("~/Downloads/PhD/Analysis/Tianeptine/Data/Covariates/IV_cov_within_A.csv")

# Append the covariate df to final_df to create the ultimate df
ult_df <- cbind(final_df, df_iv_cov)

# Fix the variables types
ult_df[["drug"]] <- as.factor(ult_df[["drug"]])
ult_df[["id"]] <- as.factor(ult_df[["id"]])

###########################

# Create a vector of target voxel names
target_voxels <- paste0("V", 1:196629)

###########################

# Define a function to fit the LMM
fit_lmm <- function(target_voxel, index) {
  # Check if all values in the target voxel column are zeros
  if(all(ult_df[[target_voxel]] == 0)) {
    # If all values are zeros, skip the column
    message(paste("Skipping", target_voxel, "because all values are zeros."))
    return(NULL)  # Return NULL to indicate skipping
  } else {
    # If not all values are zeros, fit the LMM model
    message(paste("Running a LMM on", target_voxel, "..."))
    model <- lmerTest::lmer(
      formula = as.formula(paste(target_voxel, "~ drug + mFD + (1|id)")),
      data = ult_df
    )
    return(list(index = index, model = model))  # Return index and fitted model
  }
}

# Use purrr::map to apply the function to each target voxel
lmm_results <- purrr::map(seq_along(target_voxels), ~fit_lmm(target_voxels[.x], .x))

###########################

# Define a function to extract results or return default values
extract_results_or_default <- function(model_results, term, index) {
  if (is.null(model_results[[index]])) {
    return(c(statistic = 0, p_value = 0, dof = 0))
  } else {
    stat_value <- summary(model_results[[index]]$model)$coef[term, "t value"]
    p_value <- summary(model_results[[index]]$model)$coef[term, "Pr(>|t|)"]
    dof <- summary(model_results[[index]]$model, ddf = "Satterthwaite")$coef[term, "df"]
    return(c(statistic = stat_value, p_value = p_value, dof = dof))
  }
}

# For drug_results
drug_results <- purrr::map(seq_along(lmm_results), function(i) {
  print(paste("Processing index:", i))
  extract_results_or_default(
    lmm_results, 
    "drug1", 
    i
  )
})

###########################################

# Create data frames for the statistics and p-values
drug_df <- data.frame(do.call(rbind, drug_results))

# Set the column names
#colnames(drug_df) <- c("statistic", "p_value")

###########

write.csv(drug_df, "drug_results_within_A.csv", row.names = FALSE)

###########

drug_df <- read.csv("~/Downloads/PhD/Analysis/Tianeptine/Static_FC/Jupyter_test_notebook/Analysis/2nd-level_GLM/Green_and_orange_dataset/wDC/drug_results_within_A.csv")

#############################

# Define path to the NIfTI file and output file names
target_nifti_file <- subj_list[1]
output_file <- "drug_statistics_within_A"

convert_and_save_nifti <- function(target_nifti_file, output_file, drug_df) {
  # Read the target NIfTI file
  target_image <- readnii(target_nifti_file)
  
  # Convert the first column of the dataframe to NIfTI array
  drug_nifti <- niftiarr(target_image, drug_df$statistic)
  # Alternatively, https://neuroconductor.org/help/neurobase/reference/remake_img.html
  
  # Save the NIfTI file
  writenii(drug_nifti, output_file)
}

# Usage
convert_and_save_nifti(target_nifti_file, output_file, drug_df)

###########################
# pTFCE without permutations, with multiple comparisons correction - https://github.com/spisakt/pTFCE/wiki/3.-R-package and https://rdrr.io/github/spisakt/pTFCE/man/ptfce.html 

setwd('/Users/mishodimitrov/Downloads/PhD/Analysis/Tianeptine/Static_FC/Jupyter_test_notebook/Analysis/2nd-level_GLM/Green_and_orange_dataset/wDC')

# Load mask
MASK=readNIfTI("corrected_int_mask.nii.gz")

###################################
# Load t-map
Tmap=readNIfTI("drug_statistics_within_A.nii.gz")

# Convert to z-map
degrees_of_freedom <- median(drug_df[drug_df$dof != 0, ]$dof)
Z=qnorm(pt(Tmap, df=degrees_of_freedom, log.p = T), log.p = T )

# Run pTFCE
pTFCE=ptfce(Z, MASK)

# Visualise original andpTFCE Z-maps
orthographic(Z, zlim=c(0, max(pTFCE$Z)), crosshair=F) #original
orthographic(pTFCE$Z, zlim=c(0, max(pTFCE$Z)), crosshair=F) #pTFCE
# ...and save them
writeNIfTI(Z, "/Users/mishodimitrov/Downloads/PhD/Analysis/Tianeptine/Static_FC/Jupyter_test_notebook/Analysis/2nd-level_GLM/Green_and_orange_dataset/wDC/drug_within_A_original-z-score-map")
writeNIfTI(pTFCE$Z, "/Users/mishodimitrov/Downloads/PhD/Analysis/Tianeptine/Static_FC/Jupyter_test_notebook/Analysis/2nd-level_GLM/Green_and_orange_dataset/wDC/drug_within_A-pTFCE-z-score-map")

# Visualise and save pTFCE p-values
orthographic(pTFCE$p, zlim=c(0, max(pTFCE$p)), crosshair=F)
# ...and save them
writeNIfTI(pTFCE$p, "/Users/mishodimitrov/Downloads/PhD/Analysis/Tianeptine/Static_FC/Jupyter_test_notebook/Analysis/2nd-level_GLM/Green_and_orange_dataset/wDC/drug_within_A_pTFCE-p-value-map")
