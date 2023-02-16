library(tidyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(car)
library(effsize)
library(grid)

### Code from:
### http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

### Drug effect

## Load and arrange the data
d <- read.csv("ASC_Pla_CBD_wDC_Z_hypo_no_missing.csv", header = TRUE)
colnames(d) <- c("Baseline", "Drug")
d$id <- 1:nrow(d)
d2 <- gather(data = d, key = "Condition", value = "wDC", Baseline, Drug)
colnames(d2) <- c("id", "Condition", "wDC")

## Check for normality
# http://www.sthda.com/english/wiki/normality-test-in-r

# Baseline
ggdensity(data = d$Baseline, 
          title = "Density plot of baseline wDC",
          xlab = "wDC")
ggqqplot(d$Baseline)
shapiro.test(d$Baseline)

# Drug
ggdensity(data = d$Drug, 
          title = "Density plot of drug wDC",
          xlab = "wDC")
ggqqplot(d$Drug)
shapiro.test(d$Drug)

## Check for equality of variances

# If n of samples = 2 and are normally distributed, use an F-test
var.test(d$Baseline, d$Drug)

# If NOT normally distributed (for n >=2)
leveneTest(wDC ~ Condition, data = d2)

## Check for outliers
# Baseline
boxplot(d$Baseline,
        ylab = "wDC")
# Check which subjects might be deemed outliers
out_base <- boxplot.stats(d$Baseline)$out
out_ind_base <- which(d$Baseline %in% c(out_base))
out_ind_base

# Drug
boxplot(d$Drug,
        ylab = "wDC")
# Check which subjects might be deemed outliers
out_drug <- boxplot.stats(d$Drug)$out
out_ind_drug <- which(d$Drug %in% c(out_drug))
out_ind_drug

## Decide whether or not to remove outliers
# Remove if necessary..

## Choose whether to use a parametric or a non-parametric test
# Run the test

compare_means(wDC ~ Condition, data = d2,
              #method = "t.test",
              method = "wilcox.test", 
              paired = TRUE)

## Calculate the effect size (https://www.statisticshowto.com/hedges-g/)
# For small samples (n < 20), it's preferable to use Hedges's G instead of Cohen's D
# If standard deviations are significantly different between groups, choose Glass’s delta instead. Glass’s delta uses only the control group’s standard deviation (SDC).

# https://www.rdocumentation.org/packages/effsize/versions/0.8.1/topics/cohen.d
effsize <- cohen.d(d$Drug, d$Baseline, hedges.correction = TRUE
                   #, pooled = FALSE
)
effsize_rounded <- round(effsize[["estimate"]], digits = 2)

# Add the effect size value
base_effsize_text <- "Effect Size ="
effsize_text <- paste(base_effsize_text, effsize_rounded)
grob <- grobTree(textGrob(effsize_text, x=0.25,  y=0.93, hjust=0,
                          gp=gpar(col="black", fontsize=10, fontface="italic")))


# Plot the results
p <- ggboxplot(d2, x = "Condition", y = "wDC",
               #color = "supp", 
               fill = "Condition",
               palette = c("#30fc03", "#c230b3"),
               add = "jitter"#,
               #ylim = c(-0.00000002, 0.0000ß0002)
               )

#  Add p-value and effect size
p + stat_compare_means(
  method = "wilcox.test"
  #method = "t.test"
  ,paired = TRUE
) + 
  ggtitle("Drug Effect on Mean Sensorimotor Weighted Degree Centrality in ASC") + 
  annotation_custom(grob)

