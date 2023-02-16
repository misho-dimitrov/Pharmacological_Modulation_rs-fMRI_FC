library(tidyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(car)
library(effsize)
library(grid)

### Code from:
### http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

### Baseline differences

## Load and arrange the data
d <- read.csv("TD_ASC_shift_wDC_Z_hypo.csv", header = TRUE)
d$id <- 1:nrow(d)
d2 <- gather(data = d, key = "Group", value = "wDC", TD, ASC)
colnames(d2) <- c("id", "Group", "wDC")

## Check for normality
# http://www.sthda.com/english/wiki/normality-test-in-r

# TD
ggdensity(data = d$TD, 
          title = "Density plot of TD wDC",
          xlab = "wDC")
ggqqplot(d$TD)
shapiro.test(d$TD)

# ASC
ggdensity(data = d$ASC, 
          title = "Density plot of ASC wDC",
          xlab = "wDC")
ggqqplot(d$ASC)
shapiro.test(d$ASC)

## Check for equality of variances


# If n of samples = 2 and are normally distributed, use an F-test
var.test(d$TD, d$ASC)

# If NOT normally distributed (for n >=2)
leveneTest(wDC ~ Group, data = d2)

## Check for outliers
# TD
boxplot(d$TD,
        ylab = "wDC")
# Check which subjects might be deemed outliers
out_td <- boxplot.stats(d$TD)$out
out_ind_td <- which(d$TD %in% c(out_td))
out_ind_td

# ASC
boxplot(d$ASC,
        ylab = "wDC")
# Check which subjects might be deemed outliers
out_asc <- boxplot.stats(d$ASC)$out
out_ind_asc <- which(d$ASC %in% c(out_asc))
out_ind_asc

## Decide whether or not to remove outliers
# Remove if necessary..

## Choose whether to use a parametric or a non-parametric test
# Run the test
compare_means(wDC ~ Group, data = d2,
              #method = "t.test",
              method = "wilcox.test",
              paired = FALSE)

## Calculate the effect size (https://www.statisticshowto.com/hedges-g/)
# For small samples (n < 20), it's preferable to use Hedges's G instead of Cohen's D
# If standard deviations are significantly different between groups, choose Glass’s delta instead. Glass’s delta uses only the control group’s standard deviation (SDC).

# https://www.rdocumentation.org/packages/effsize/versions/0.8.1/topics/cohen.d
effsize <- cohen.d(d$ASC, d$TD, hedges.correction = TRUE
                   #, pooled = FALSE
                   )
effsize_rounded <- round(effsize[["estimate"]], digits = 2)

## Plot the results
p <- ggboxplot(d2, x = "Group", y = "wDC",
               #color = "supp", 
               fill = "Group",
               palette = c("#E7B800", "#00AFBB"),
               add = "jitter") 

# Add the effect size value
base_effsize_text <- "Effect Size ="
effsize_text <- paste(base_effsize_text, effsize_rounded)
grob <- grobTree(textGrob(effsize_text, x=0.245,  y=0.93, hjust=0,
                          gp=gpar(col="black", fontsize=10, fontface="italic")))

#  Add p-value and effect size
p + stat_compare_means(
  #method = "t.test"
  method = "wilcox.test"
  #, label.x = 0.7
  #, label.y = 0.77
  ) + 
  ggtitle("Shift Differences in Mean Sensorimotor Weighted Degree Centrality") + 
  annotation_custom(grob)

