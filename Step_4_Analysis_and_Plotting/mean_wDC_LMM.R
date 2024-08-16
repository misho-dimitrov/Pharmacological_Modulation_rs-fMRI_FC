# Load required libraries
library(ptestr)
library(readr)

setwd('~/Downloads/PhD/Analysis/Tianeptine/Static_FC/Jupyter_test_notebook/Analysis/Mean_wDC/')

########################################################

### Comprehensive models ###
hyper_wDC_df <- read.csv("hyper_LMM.csv")
hyper_wDC_df <- hyper_wDC_df[, -1]

hyper_wDC.permuted <- grouped_perm_glmm(hyper_wDC_df, formla = wDC ~ group * drug + mFD + (1|ID), var_to_perm='wDC', permNum = 5000, seed = 42)

write_csv(hyper_wDC.permuted, "hyper_model_output.csv")
###############
hypo_wDC_df <- read.csv("hypo_LMM.csv")
hypo_wDC_df <- hypo_wDC_df[, -1]

hypo_wDC.permuted <- grouped_perm_glmm(hypo_wDC_df, formla = wDC ~ group * drug + mFD + (1|ID), var_to_perm='wDC', permNum = 5000, seed = 42)

write_csv(hypo_wDC.permuted, "hypo_model_output.csv")

########################################################

### Within-group Models ###
# Filter rows where group = 0
hyper_wDC_df_NT <- subset(hyper_wDC_df, group == 0)

# Remove the group column
hyper_wDC_df_NT <- hyper_wDC_df_NT[, -which(names(hyper_wDC_df_NT) == "group")]

# Alternatively, you can remove the column by column index
# hyper_wDC_df_NT <- hyper_wDC_df_NT[, -2]

hyper_wDC_NT.permuted <- grouped_perm_glmm(hyper_wDC_df_NT, formla = wDC ~ drug + mFD + (1|ID), var_to_perm='wDC', permNum = 5000, seed = 42)

write_csv(hyper_wDC_NT.permuted, "hyper_NT_model_output.csv")
###############
# Filter rows where group = 0
hyper_wDC_df_A <- subset(hyper_wDC_df, group == 1)

# Remove the group column
hyper_wDC_df_A <- hyper_wDC_df_A[, -which(names(hyper_wDC_df_A) == "group")]

# Alternatively, you can remove the column by column index
# hyper_wDC_df_A <- hyper_wDC_df_A[, -2]

hyper_wDC_A.permuted <- grouped_perm_glmm(hyper_wDC_df_A, formla = wDC ~ drug + mFD + (1|ID), var_to_perm='wDC', permNum = 5000, seed = 42)

write_csv(hyper_wDC_A.permuted, "hyper_A_model_output.csv")
###############
# Filter rows where group = 0
hypo_wDC_df_NT <- subset(hypo_wDC_df, group == 0)

# Remove the group column
hypo_wDC_df_NT <- hypo_wDC_df_NT[, -which(names(hypo_wDC_df_NT) == "group")]

# Alternatively, you can remove the column by column index
# hypo_wDC_df_NT <- hypo_wDC_df_NT[, -2]

hypo_wDC_NT.permuted <- grouped_perm_glmm(hypo_wDC_df_NT, formla = wDC ~ drug + mFD + (1|ID), var_to_perm='wDC', permNum = 5000, seed = 42)

write_csv(hypo_wDC_NT.permuted, "hypo_NT_model_output.csv")
###############
# Filter rows where group = 0
hypo_wDC_df_A <- subset(hypo_wDC_df, group == 1)

# Remove the group column
hypo_wDC_df_A <- hypo_wDC_df_A[, -which(names(hypo_wDC_df_A) == "group")]

# Alternatively, you can remove the column by column index
# hypo_wDC_df_A <- hypo_wDC_df_A[, -2]

hypo_wDC_A.permuted <- grouped_perm_glmm(hypo_wDC_df_A, formla = wDC ~ drug + mFD + (1|ID), var_to_perm='wDC', permNum = 5000, seed = 42)

write_csv(hypo_wDC_A.permuted, "hypo_A_model_output.csv")
