library(tidyr)
library(ggplot2)
library(ggpubr)

### http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r

d <- read.csv("5-HT_levels_AQ.csv", header = TRUE#, colClasses=c(NA,"NULL","NULL", NA)
              )
colnames(d) <- c("Serotonin", "AQ")

d_prelim <- ggplot(d, aes(x=Serotonin, y=AQ#, col=Group
                          )) + geom_point()
annotate_figure(d_prelim, top = text_grob("Scatterplot of Blood Serotonin Levels and AQ", 
                                         #color = "red", face = "bold", size = 14
))

## Check for normality
# http://www.sthda.com/english/wiki/normality-test-in-r

# AQ
ggdensity(data = d$AQ, 
          title = "Density plot of shift AQ",
          xlab = "AQ")
ggqqplot(d$AQ)
shapiro.test(d$AQ)

# 5-HT
ggdensity(data = d$Serotonin, 
          title = "Density plot of Serotonin",
          xlab = "Serotonin")
ggqqplot(d$Serotonin)
shapiro.test(d$Serotonin)

## Do the correlation analysis and plot the results
sp <- ggscatter(d, x = "Serotonin", y = "AQ",
                add = "reg.line",  # Add regression line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                xlab = "Blood Serotonin (ng/mL)") + 
  geom_vline(xintercept=270, linetype="dashed", color = "red")

sp

# Add correlation coefficient
sp_corr <- sp +
  stat_cor(#method = "pearson" 
              method = "spearman"
              #,label.x = 3, label.y = 30
              )


annotate_figure(sp_corr, top = text_grob("Correlation between Blood Serotonin Levels and AQ", 
                                              #color = "red", face = "bold", size = 14
))