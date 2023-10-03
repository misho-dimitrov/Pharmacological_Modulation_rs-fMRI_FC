#################################                  #####################################
############### !!!!!!!!!!!!!!!!! CREATE FUNCTIONS !!!!!!!!!!!!!!!!! ###################
                # (too much (necessary and unnecessary) repetition..)
#################################                  #####################################


### Load the libraries
library(tidyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(car)
library(stringr)
library(effsize)
library(grid)
library(yhat)
library(gvlma)
library(lmerTest)
library(lattice)
library(data.table)

setwd('~/Downloads/PhD/Analysis/Mean_wDC/Data/')

## Create a list of csv files
csv_list <- c("ARB_baseline_hyper.csv", "ARB_baseline_hypo.csv", "ARB_drug_effect_ASC_hyper.csv", "ARB_drug_effect_ASC_hypo.csv", "ARB_drug_effect_TD_hyper.csv", "ARB_drug_effect_TD_hypo.csv", "CBD_baseline_hyper.csv", "CBD_baseline_hypo.csv", "CBD_drug_effect_ASC_hyper.csv", "CBD_drug_effect_ASC_hypo.csv", "CBD_drug_effect_TD_hyper.csv", "CBD_drug_effect_TD_hypo.csv", "CIT_baseline_hyper.csv", "CIT_baseline_hypo.csv", "CIT_drug_effect_ASC_hyper.csv", "CIT_drug_effect_ASC_hypo.csv", "CIT_drug_effect_TD_hyper.csv", "CIT_drug_effect_TD_hypo.csv", "TIA_baseline_hyper.csv", "TIA_baseline_hypo.csv", "TIA_drug_effect_ASC_hyper.csv", "TIA_drug_effect_ASC_hypo.csv", "TIA_drug_effect_TD_hyper.csv", "TIA_drug_effect_TD_hypo.csv")
stats_csv_list <- c()

p_val_df <- data.frame(matrix(ncol = 2, nrow = 1))
p_val_df_cols <- c("comparison_name", "p_value")
colnames(p_val_df) <- p_val_df_cols

# Iterate over each csv file
for (csv in csv_list) {
  # Load the data
  d <- read.csv(csv, header = TRUE#, na.rm = TRUE
                )
  # Create work variables
  paired <- d$Paired[1]
  within_g <- d$Within[1]
  mask <- d$Mask[1]
  balanced <- d$Balanced[1]
  
  # Use a t-test/Wicoxon test if balanced or a linear effects model with mFD as a covariate if imbalanced
  if(balanced == TRUE) {
    d2 <- select(d, wDC, Group)
    id_group1 <- 1:nrow(d2[d2$Group == "0",])
    id_group2 <- 1:nrow(d2[d2$Group == "1",])
    d2$id <- c(id_group1, id_group2)
    
    # Check if paired or unpaired before testing
    if(paired == FALSE) {
      colnames(d2) <- c("wDC", "Group", "id")
      d2$Group[d2$Group == "0"] <- "Non-autistic"
      d2$Group[d2$Group == "1"] <- "ASC"
      
      ## Check for normality
      # http://www.sthda.com/english/wiki/normality-test-in-r
      # Non-autistic
      ggdensity(data = d2$wDC[d2$Group == "Non-autistic"], 
                title = "Density plot of Non-autistic",
                xlab = "wDC")
      ggqqplot(d2$wDC[d2$Group == "Non-autistic"])
      shapiro_TD_p <- shapiro.test(d2$wDC[d2$Group == "Non-autistic"])$p.value
      # ASC
      ggdensity(data = d2$wDC[d2$Group == "ASC"], 
                title = "Density plot of ASC",
                xlab = "wDC")
      ggqqplot(d2$wDC[d2$Group == "ASC"])
      shapiro_ASC_p <- shapiro.test(d2$wDC[d2$Group == "ASC"])$p.value
      
      if((shapiro_TD_p > 0.05) & (shapiro_ASC_p > 0.05)) {
        norm_dist <- TRUE
      } else {
        norm_dist <- FALSE
      }
      
      ## Check for equality of variances
      if(norm_dist == TRUE){
        bartlett_test <- bartlett.test(wDC ~ Group, data = d2)
        
        if(bartlett_test$p.value > 0.05) {
          eq_var = TRUE 
        } else {
          eq_var = FALSE
        }
        
      } else if(norm_dist == FALSE) {
        levene_test <- leveneTest(wDC ~ Group, data = d2)
        
        if(levene_test$`Pr(>F)` > 0.05) {
          eq_var = TRUE 
        } else {
          eq_var = FALSE
        }
        
      }
      
      ## Get summary statistics
      mean_TD <- mean(d2$wDC[d2$Group == "Non-autistic"])
      mean_ASC <- mean(d2$wDC[d2$Group == "ASC"])
      median_TD <- median(d2$wDC[d2$Group == "Non-autistic"])
      median_ASC <- median(d2$wDC[d2$Group == "ASC"])
      SD_TD <- sd(d2$wDC[d2$Group == "Non-autistic"])
      SD_ASC <- sd(d2$wDC[d2$Group == "ASC"])
      IQR_TD <- IQR(d2$wDC[d2$Group == "Non-autistic"])
      IQR_ASC <- IQR(d2$wDC[d2$Group == "ASC"])

      ## Run statistical analysis
      if(norm_dist == TRUE) {
        test_method <- "t.test"
      } else if(norm_dist == FALSE) {
        test_method <- "wilcox.test"
      }
    
      stat_test <- compare_means(wDC ~ Group, data = d2,
                                 method = test_method, 
                                 paired = paired)
      
      stat_test_method <- stat_test$method
      p_val <- stat_test$p
      
      ## Calculate the effect size (https://www.statisticshowto.com/hedges-g/)
      # For small samples (n < 20), it's preferable to use Hedges's G instead of Cohen's D
      # If standard deviations are significantly different between groups, choose Glass’s delta instead. Glass’s delta uses only the control group’s standard deviation (SDC).
      
      # https://www.rdocumentation.org/packages/effsize/versions/0.8.1/topics/cohen.d
      if(eq_var == TRUE ) {
        effsize <- cohen.d(d2$wDC[d2$Group == "ASC"], d2$wDC[d2$Group == "Non-autistic"], hedges.correction = TRUE
                           #, pooled = FALSE
        )
        effsize_rounded <- round(effsize[["estimate"]], digits = 2)
      } else if(eq_var == FALSE) {
        effsize <- cohen.d(d2$wDC[d2$Group == "ASC"], d2$wDC[d2$Group == "Non-autistic"], hedges.correction = TRUE
                           , pooled = FALSE
        )
        effsize_rounded <- round(effsize[["estimate"]], digits = 2)
      }
      
      ## Save all of the statistical information in a list and on disk
      stat_info <- rbind(csv, 
                         paired, 
                         within_g, 
                         mask, 
                         balanced, 
                         shapiro_TD_p,
                         shapiro_ASC_p,
                         norm_dist, 
                         eq_var, 
                         mean_TD,
                         mean_ASC,
                         median_TD,
                         median_ASC,
                         SD_TD,
                         SD_ASC,
                         IQR_TD,
                         IQR_ASC,
                         stat_test_method,
                         p_val,
                         effsize_rounded)
      
      stats_csv_pre <- str_split(csv, "\\.")[[1]][1]
      stats_csv <- paste(stats_csv_pre,"_stats.csv", sep = "")
      stats_csv_list <- c(stats_csv_list, stats_csv)
      stats_csv_dir <- paste("~/Downloads/PhD/Analysis/Mean_wDC/Data/", stats_csv, sep = "")
      write.csv(stat_info, stats_csv_dir, row.names=TRUE)
      
      # Add the p-value to the p-value dataframe
      p_val_df <- rbind(p_val_df, c(stats_csv_pre, p_val))
      
    } else if(paired == TRUE) {
      colnames(d2) <- c("wDC", "Condition", "id")
      d2$Condition[d2$Condition == "0"] <- "Baseline"
      d2$Condition[d2$Condition == "1"] <- "Drug"
      
      ## Check for normality
      # http://www.sthda.com/english/wiki/normality-test-in-r
      # Baseline
      ggdensity(data = d2$wDC[d2$Condition == "Baseline"], 
                title = "Density plot of Baseline",
                xlab = "wDC")
      ggqqplot(d2$wDC[d2$Condition == "Baseline"])
      shapiro_TD_p <- shapiro.test(d2$wDC[d2$Condition == "Baseline"])$p.value
      # Drug
      ggdensity(data = d2$wDC[d2$Condition == "Drug"], 
                title = "Density plot of Drug",
                xlab = "wDC")
      ggqqplot(d2$wDC[d2$Condition == "Drug"])
      shapiro_ASC_p <- shapiro.test(d2$wDC[d2$Condition == "Drug"])$p.value
      
      if((shapiro_TD_p > 0.05) & (shapiro_ASC_p > 0.05)) {
        norm_dist <- TRUE
      } else {
        norm_dist <- FALSE
      }
      
      ## Check for equality of variances
      if(norm_dist == TRUE){
        bartlett_test <- bartlett.test(wDC ~ Condition, data = d2)
        
        if(bartlett_test$p.value > 0.05) {
          eq_var = TRUE 
        } else {
          eq_var = FALSE
        }
        
      } else if(norm_dist == FALSE) {
        levene_test <- leveneTest(wDC ~ Condition, data = d2)
        
        if(levene_test$`Pr(>F)` > 0.05) {
          eq_var = TRUE 
        } else {
          eq_var = FALSE
        }
        
      }
      
      ## Get summary statistics
      mean_Baseline <- mean(d2$wDC[d2$Condition == "Baseline"])
      mean_Drug <- mean(d2$wDC[d2$Condition == "Drug"])
      median_Baseline <- median(d2$wDC[d2$Condition == "Baseline"])
      median_Drug <- median(d2$wDC[d2$Condition == "Drug"])
      SD_Baseline <- sd(d2$wDC[d2$Condition == "Baseline"])
      SD_Drug <- sd(d2$wDC[d2$Condition == "Drug"])
      IQR_Baseline <- IQR(d2$wDC[d2$Condition == "Baseline"])
      IQR_Drug <- IQR(d2$wDC[d2$Condition == "Drug"])
      
      ## Run statistical analysis
      if(norm_dist == TRUE) {
        test_method <- "t.test"
      } else if(norm_dist == FALSE) {
        test_method <- "wilcox.test"
      }
      
      stat_test <- compare_means(wDC ~ Condition, data = d2,
                                 method = test_method, 
                                 paired = paired)
      
      stat_test_method <- stat_test$method
      p_val <- stat_test$p
      
      ## Calculate the effect size (https://www.statisticshowto.com/hedges-g/)
      # For small samples (n < 20), it's preferable to use Hedges's G instead of Cohen's D
      # If standard deviations are significantly different between groups, choose Glass’s delta instead. Glass’s delta uses only the control group’s standard deviation (SDC).
      
      # https://www.rdocumentation.org/packages/effsize/versions/0.8.1/topics/cohen.d
      if(eq_var == TRUE ) {
        effsize <- cohen.d(d2$wDC[d2$Condition == "Drug"], d2$wDC[d2$Condition == "Baseline"], hedges.correction = TRUE
                           #, pooled = FALSE
        )
        effsize_rounded <- round(effsize[["estimate"]], digits = 2)
      } else if(eq_var == FALSE) {
        effsize <- cohen.d(d2$wDC[d2$Condition == "Drug"], d2$wDC[d2$Condition == "Baseline"], hedges.correction = TRUE
                           , pooled = FALSE
        )
        effsize_rounded <- round(effsize[["estimate"]], digits = 2)
      }
      
      ## Save all of the statistical information in a list and on disk
      stat_info <- rbind(csv, 
                         paired, 
                         within_g, 
                         mask, 
                         balanced, 
                         norm_dist, 
                         eq_var, 
                         mean_Baseline,
                         mean_Drug,
                         median_Baseline,
                         median_Drug,
                         SD_Baseline,
                         SD_Drug,
                         IQR_Baseline,
                         IQR_Drug,
                         stat_test_method,
                         p_val,
                         effsize_rounded)
      
      stats_csv_pre <- str_split(csv, "\\.")[[1]][1]
      stats_csv <- paste(stats_csv_pre,"_stats.csv", sep = "")
      stats_csv_list <- c(stats_csv_list, stats_csv)
      stats_csv_dir <- paste("~/Downloads/PhD/Analysis/Mean_wDC/Data/", stats_csv, sep = "")
      write.csv(stat_info, stats_csv_dir, row.names=TRUE)
      
      # Add the p-value to the p-value dataframe
      p_val_df <- rbind(p_val_df, c(stats_csv_pre, p_val))
      
    } 
    
  } else if(balanced == FALSE) {
    d2 <- select(d, wDC, Group, mFD)
    id_group1 <- 1:nrow(d2[d2$Group == "0",])
    id_group2 <- 1:nrow(d2[d2$Group == "1",])
    d2$id <- c(id_group1, id_group2)
    
    if(paired == FALSE) {
      colnames(d2) <- c("wDC", "Group", "mFD", "id")
      d2$Group[d2$Group == "0"] <- "Non-autistic"
      d2$Group[d2$Group == "1"] <- "ASC"
      
      ## Get summary statistics
      mean_TD <- mean(d2$wDC[d2$Group == "Non-autistic"])
      mean_ASC <- mean(d2$wDC[d2$Group == "ASC"])
      median_TD <- median(d2$wDC[d2$Group == "Non-autistic"])
      median_ASC <- median(d2$wDC[d2$Group == "ASC"])
      SD_TD <- sd(d2$wDC[d2$Group == "Non-autistic"])
      SD_ASC <- sd(d2$wDC[d2$Group == "ASC"])
      IQR_TD <- IQR(d2$wDC[d2$Group == "Non-autistic"])
      IQR_ASC <- IQR(d2$wDC[d2$Group == "ASC"])
      
      ## Perform the linear regression with mFD as a covariate
      # Set up the variables
      stat_test_method <- "mFD-adjusted T-test"
      wDC <- unlist(d['wDC'])
      Group <- unlist(d['Group'])
      mFD <- unlist(d['mFD'])
      
      # If you run the following model, you will get the same inference as if you do a t-test
      model_simple <- lm(wDC ~ Group)
      summary(model_simple)
      # Run the model with mFD as covariate
      model_with_cov <- lm(wDC ~ Group + mFD)
      summary(model_with_cov)
      
      # Check whether the model assumptions are met
      gvlma_object <- gvlma(model_with_cov)
      lm_assumptions <- gvlma_object$GlobalTest$GlobalStat4$Decision
      
      if(gvlma_object$GlobalTest$GlobalStat4$Decision == 0) {
        # Define function to extract overall p-value of model
        overall_p <- function(my_model) {
          f <- summary(my_model)$fstatistic
          p <- pf(f[1],f[2],f[3],lower.tail=F)
          attributes(p) <- NULL
          p_rounded <- round(p, digits = 2)
          return(p_rounded)
        }
        
        # Extract overall p-value of model
        p_val <- overall_p(model_with_cov)
        
        # Calculate effect size using different algorithms and save the recommended one
        effsize <- effect.size(model_with_cov)
        effsize_recommended <- as.double(effsize[effsize$Recommended == 'Yes',]$Effect.Size)
        effsize_rounded <- round(effsize_recommended, digits = 2)
        
      } else if(gvlma_object$GlobalTest$GlobalStat4$Decision != 0){
        print("Regression assumptions violated. Use at your own risk..")
        
        # Define function to extract overall p-value of model
        overall_p <- function(my_model) {
          f <- summary(my_model)$fstatistic
          p <- pf(f[1],f[2],f[3],lower.tail=F)
          attributes(p) <- NULL
          p_rounded <- round(p, digits = 2)
          return(p_rounded)
        }
        
        # Extract overall p-value of model
        p_val <- overall_p(model_with_cov)
        
        # Calculate effect size using different algorithms and save the recommended one
        effsize <- effect.size(model_with_cov)
        effsize_recommended <- as.double(effsize[effsize$Recommended == 'Yes',]$Effect.Size)
        effsize_rounded <- round(effsize_recommended, digits = 2)
      }
      
      ## Save all of the statistical information in a list and on disk
      stat_info <- rbind(csv, 
                         paired, 
                         within_g, 
                         mask, 
                         balanced, 
                         mean_TD,
                         mean_ASC,
                         median_TD,
                         median_ASC,
                         SD_TD,
                         SD_ASC,
                         IQR_TD,
                         IQR_ASC,
                         stat_test_method,
                         lm_assumptions,
                         p_val,
                         effsize_rounded)
      
      stats_csv_pre <- str_split(csv, "\\.")[[1]][1]
      stats_csv <- paste(stats_csv_pre,"_stats.csv", sep = "")
      stats_csv_list <- c(stats_csv_list, stats_csv)
      stats_csv_dir <- paste("~/Downloads/PhD/Analysis/Mean_wDC/Data/", stats_csv, sep = "")
      write.csv(stat_info, stats_csv_dir, row.names=TRUE)
      
      # Add the p-value to the p-value dataframe
      p_val_df <- rbind(p_val_df, c(stats_csv_pre, p_val))
      
    } else if(paired == TRUE) { 
      colnames(d2) <- c("wDC", "Condition", "mFD", "id")
      d2$Condition[d2$Condition == "0"] <- "Baseline"
      d2$Condition[d2$Condition == "1"] <- "Drug"
      
      ## Get summary statistics
      mean_Baseline <- mean(d2$wDC[d2$Condition == "Baseline"])
      mean_Drug <- mean(d2$wDC[d2$Condition == "Drug"])
      median_Baseline <- median(d2$wDC[d2$Condition == "Baseline"])
      median_Drug <- median(d2$wDC[d2$Condition == "Drug"])
      SD_Baseline <- sd(d2$wDC[d2$Condition == "Baseline"])
      SD_Drug <- sd(d2$wDC[d2$Condition == "Drug"])
      IQR_Baseline <- IQR(d2$wDC[d2$Condition == "Baseline"])
      IQR_Drug <- IQR(d2$wDC[d2$Condition == "Drug"])
      
      ## Perform the LME with mFD as a covariate
      # Set up the variables
      stat_test_method <- "mFD-adjusted T-test"
      #wDC <- unlist(d['wDC'])
      #Condition <- unlist(d['Condition'])
      #mFD <- unlist(d['mFD'])
      
      # Run the LME model with mFD as covariate
      # https://ademos.people.uic.edu/Chapter18.html
      # https://benwhalley.github.io/just-enough-r/extending-traditional-rm-anova.html
      wDC.lmer <- lmer(wDC ~ Condition + mFD + (1|id), data=d2, na.action=na.exclude)
      
      # Check the summary results of the model
      summary(wDC.lmer)
      model_summary <- anova(wDC.lmer)
      
      ## Check the model assumptions
      # Linearity
      Plot.Model.wDC.lmer.Linearity <- plot(resid(wDC.lmer),d2$wDC)
      #orig <- d2$wDC[!is.na(d2$wDC)]
      #res <- resid(wDC.lmer)
      #orig_vs_res <- cbind(orig, res)
      #linearity_df <- data.table::as.data.table(orig_vs_res)
      #ggplot(linearity_df, aes(x=orig, y=res)) + geom_point()
      
      # Homogeneity of variance
      wDC.lmer.Res <- residuals(wDC.lmer) #extracts the residuals and places them in a new column in our original data table
      Abs.wDC.lmer.Res <-abs(wDC.lmer.Res) #creates a new column with the absolute value of the residuals
      wDC.lmer.Res2 <- wDC.lmer.Res^2 #squares the absolute values of the residuals to provide the more robust estimate
      wDC.lmer.Res2_df <- data.table::as.data.table(wDC.lmer.Res2, keep.rownames=TRUE)
      Levene.Model.F <- lm(wDC.lmer.Res2 ~ rn, data=wDC.lmer.Res2_df) #ANOVA of the squared residuals
      anova(Levene.Model.F) #displays the results
      
      # Normal distribution of residuals
      qqmath(wDC.lmer, id=0.05)
      
      # Get the p-value
      p_val <- model_summary['Condition', 'Pr(>F)']
      
      # Calculate effect size .. or not?
      # See https://stats.stackexchange.com/questions/546957/is-there-any-function-in-r-that-computes-effect-sizes-based-on-linear-mixed-effe
      # The R emmeans package provides a way to calculate Cohen's d from many types of models when it makes sense, with an eff_size() function. See the pairwise comparisons section of the vignette for an introduction.
      # Nevertheless, the package author, Russ Lenth, discusses why presenting Cohen's d isn't always a good idea, in comments on this answer about the degrees of freedom (d.f) to specify to get an "effect size" for a mixed model. In particular:
      # In a simple lm() model, things are a lot more straightforward. In a mixed model, there is a lot more ambiguity. But to me, the d.f. to use is the smallest issue. What's bigger is deciding what you are even talking about when you compute an effect size. There is a lot of discussion on this, and from where I sit, the question is nearly unanswerable. I don't really believe in effect sizes. I provided the function because I thought [it] was important to account for uncertainty in the SD estimate if you insist on computing an effect size.
      effsize_rounded <- NaN
      
      ## Save all of the statistical information in a list and on disk
      stat_info <- rbind(csv, 
                         paired, 
                         within_g, 
                         mask, 
                         balanced, 
                         mean_Baseline,
                         mean_Drug,
                         median_Baseline,
                         median_Drug,
                         SD_Baseline,
                         SD_Drug,
                         IQR_Baseline,
                         IQR_Drug,
                         stat_test_method,
                         p_val,
                         effsize_rounded
                         )
      
      stats_csv_pre <- str_split(csv, "\\.")[[1]][1]
      stats_csv <- paste(stats_csv_pre,"_stats.csv", sep = "")
      stats_csv_list <- c(stats_csv_list, stats_csv)
      stats_csv_dir <- paste("~/Downloads/PhD/Analysis/Mean_wDC/Data/", stats_csv, sep = "")
      write.csv(stat_info, stats_csv_dir, row.names=TRUE)
      
      # Add the p-value to the p-value dataframe
      p_val_df <- rbind(p_val_df, c(stats_csv_pre, p_val))
    }
   
  } 
  
  # Clear the workspace except the p-value dataframe
  rm(list=setdiff(ls(), c("csv_list", "stats_csv_list", "p_val_df")))
  
}

## Get the p-value dataframe and correct the p-values for multiple comparisons
p_val_list <- p_val_df[c(2:length(p_val_df[, "p_value"])), "p_value"]
corr_p_val_list <- p.adjust(p_val_list, method = "BY")
# Store the corrected p-values in a new dataframe similar to the old one
corr_p_val_df <- data.frame(p_val_df)
corr_p_val_df[c(2:length(corr_p_val_df[, "comparison_name"])), "comparison_name"] <- stats_csv_list
corr_p_val_df[c(2:length(corr_p_val_df[, "p_value"])), "p_value"] <- corr_p_val_list
# Save to disk
#corr_p_val_df_comparison <- str_split(corr_p_val_df[2,1], "\\_")[[1]][2]
#corr_p_val_df_mask <- str_split(corr_p_val_df[2,1], "\\_")[[1]][3]
#corr_p_val_df_name <- paste("corrected_p-values_", paste(corr_p_val_df_comparison, paste("_", paste(corr_p_val_df_mask, ".csv", sep = ""), sep = ""), sep = ""), sep = "")
write.csv(corr_p_val_df, "~/Downloads/PhD/Analysis/Mean_wDC/Data/corrected_p-values.csv", row.names=TRUE)

rm(list=ls())

######                                                                            ######
########### -----------------------------------------------------------------###########
################################## PLOTTING ############################################
########### -----------------------------------------------------------------###########
######                                                                            ######

setwd('~/Downloads/PhD/Analysis/Mean_wDC/Data/')

i <- 0

### Change this manually !!!!
custom_csv_list <- c("ARB_drug_effect_ASC_hypo.csv", "CIT_drug_effect_ASC_hypo.csv", "TIA_drug_effect_ASC_hypo.csv", "CBD_drug_effect_ASC_hypo.csv")
### Change this manually !!!!
  
## Set global minimum and maximum
for (csv in custom_csv_list){
  # Set up a counter
  i <- i + 1
  
  # Load the wDC data
  d_wDC_pre <- read.csv(csv, header = TRUE)
  d_wDC <- select(d_wDC_pre, wDC, Group)
  assign(paste("d_wDC", i, sep = "_"), d_wDC)
}

# Calculate global y-axis range
global_min <- min(c(d_wDC_1$wDC, d_wDC_2$wDC, d_wDC_3$wDC, d_wDC_4$wDC))
global_max <- max(c(d_wDC_1$wDC, d_wDC_2$wDC, d_wDC_3$wDC, d_wDC_4$wDC)) + 0.2
custom_ylim <- c(global_min, global_max)

rm(list=setdiff(ls(), c("custom_csv_list", "custom_ylim")))

# Generate all subplots
i <- 0
for (csv in custom_csv_list){
  # Set up a counter
  i <- i + 1
  
  # Load the wDC data
  d_wDC_pre <- read.csv(csv, header = TRUE)
  d_wDC <- select(d_wDC_pre, wDC, Group)
  assign(paste("d_wDC", i, sep = "_"), d_wDC)
  
  # Load the stats info
  stats_csv_pre <- str_split(csv, "\\.")[[1]][1]
  stats_csv <-paste(stats_csv_pre,"_stats.csv", sep = "")
  d_stats_raw <- read.csv(stats_csv, header = TRUE, row.names=1)
  d_stats <- transpose(d_stats_raw)
  colnames(d_stats) <- rownames(d_stats_raw)
  
  # Set variables depending on the pair type
  if(d_stats$paired == FALSE) {
    title_paired <- "Baseline Differences in"
    x_label <- "Group"
    d_wDC$Group[d_wDC$Group == "0"] <- "Non-autistic"
    d_wDC$Group[d_wDC$Group == "1"] <- "Autistic"
    plot_palette <- c("#E7B800", "#00AFBB")
    stat_test_method <- d_stats$stat_test_method
    p_value <- round(as.numeric(d_stats$p_val), digits = 2)
    effsize_rounded <- d_stats$effsize_rounded
    
  } else if(d_stats$paired == TRUE) {
    title_paired <- "Drug effect on"
    title_group_pre <- str_split(d_stats$csv, "\\_")[[1]][4] # autistic vs. non-autistic
    if(title_group_pre == "TD") {
      title_group <- "Non-autistic"
    } else if(title_group_pre == "ASC") {
      title_group <- "Autistic"
    }
    colnames(d_wDC) <- c("wDC", "Condition")
    x_label <- "Condition"
    d_wDC$Condition[d_wDC$Condition == "0"] <- "Baseline"
    d_wDC$Condition[d_wDC$Condition == "1"] <- "Drug"
    plot_palette <- c("#30fc03", "#c230b3")
    stat_test_method <- d_stats$stat_test_method
    p_value <- round(as.numeric(d_stats$p_val), digits = 2)
    effsize_rounded <- d_stats$effsize_rounded
    
  } 
  
  # Set more variables depending on the mask type
  if(d_stats$mask == "Hyper") {
    title_mask <- "Fronto-parietal wDC"
    #custom_ylim <- c(-0.4, 1)
    #stat_test_label_x <- 0.645
    #stat_test_label_y <- 1
    grob1_x <- 0.25
    grob1_y <- 0.96
    grob2_x <- 0.25
    grob2_y <- 0.88
  } else if(d_stats$mask == "Hypo") {
    title_mask <- "Sensori-motor wDC"
    #custom_ylim <- c(-0.55, 1.2)
    #stat_test_label_x <- 0.645
    #stat_test_label_y <- 1.12
    grob1_x <- 0.24
    grob1_y <- 0.96
    grob2_x <- 0.25
    grob2_y <- 0.88
  } 
  
  # Choose subplot title
  if(i == 1){
    subplot_title <- "Arbaclofen"
  } else if(i == 2){
    subplot_title <- "Citalopram"
  } else if(i == 3) {
    subplot_title <- "Tianeptine"
  } else if(i == 4) {
    subplot_title <- "Cannabidiol"
  }
  
  # Create a subplot
  subplot <- ggboxplot(d_wDC, x = x_label, y = "wDC",
                         #color = "supp", 
                         fill = x_label,
                         palette = plot_palette,
                         add = "jitter") + 
    scale_y_continuous(limits=custom_ylim) + 
    ggtitle(subplot_title)
  
  # Set up the p-value bit
  #corr_p_val_csv <- read.csv("corrected_p-values.csv", header = TRUE)
  #corr_p_val <- corr_p_val_csv$p_value[corr_p_val_csv$comparison_name == stats_csv][2]
  base_p_val_text <- paste(d_stats$stat_test_method, ",")
  p_val_text <- paste(base_p_val_text, paste("p-value =", p_value)) ## USE UNADJUSTED P-VALUES
  #p_val_text <- paste(base_p_val_text, paste("p-value =", corr_p_val)) ## USE ADJUSTED P-VALUES
  grob1 <-  grobTree(textGrob(p_val_text, grob1_x,  grob1_y, hjust=0,
                              gp=gpar(col="black", fontsize=10, fontface="italic")))
  
  # Set up the effect size bit
  base_effsize_text <- "Effect Size ="
  effsize_text <- paste(base_effsize_text, effsize_rounded)
  grob2 <- grobTree(textGrob(effsize_text, grob2_x,  grob2_y, hjust=0,
                             gp=gpar(col="black", fontsize=10, fontface="italic")))
  
  #  Add p-value and effect size to plot
  final_subplot <- subplot + 
    #ggtitle(plot_title) +
    #theme(plot.title = element_text(hjust = 0.5)) + 
    annotation_custom(grob1) + 
    annotation_custom(grob2)
  
  assign(paste("subplot", i, sep = "_"), final_subplot)
}

# Add master plot title
if(d_stats$paired == FALSE){
  plot_title <- paste(title_paired, title_mask)
} else if(d_stats$paired == TRUE) {
  plot_title <- paste(title_paired, title_mask, "in", title_group, "Individuals")
}

### Combine similar plots across drugs
# https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

final_plot <- ggarrange(subplot_1, subplot_2, subplot_3, subplot_4, ncol = 2, nrow = 2, 
          common.legend = TRUE, legend="bottom") + 
  ggtitle(plot_title) +
  theme(plot.title = element_text(hjust = 0.5)) 

### SAVE EACH PLOT
ggsave(paste(plot_title, ".png", sep = ""), plot = final_plot, 
       #device = "png", dpi = 800, 
       width = 3000, height = 1800, units = "px"
       )

