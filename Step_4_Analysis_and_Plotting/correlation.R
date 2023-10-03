library(tidyr)
library(ggplot2)
library(ggpubr)

### http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r

d <- read.csv("TD_ASD_TIA_SHIFT_with_AQ_hyper.csv", header = TRUE#, colClasses=c("NULL",NA,NA)
              )

#d$Group[d$Group == "TD"] = "Non-autistic"
#d$Group[d$Group == "ASC"] = "Autistic"

## Check for linearity
ggplot(d, aes(x=wDC, y=AQ)) + geom_point()

## Check for normality
# http://www.sthda.com/english/wiki/normality-test-in-r

# wDC
ggdensity(data = d$wDC, 
          title = "Density plot of baseline wDC",
          xlab = "wDC")
ggqqplot(d$wDC)
shapiro.test(d$wDC)

# AQ
ggdensity(data = d$AQ, 
          title = "Density plot of AQ",
          xlab = "AQ")
ggqqplot(d$AQ)
shapiro.test(d$AQ)

## Check for outliers
# get from the other script

## Do the correlation analysis and plot the results
ggplot(d, aes(wDC,AQ)) + geom_point(aes(colour=Group)) + 
  scale_color_manual(values = c(`Non-autistic` = "orange", Autistic = "green")) +
  stat_smooth(method = "lm",
              geom = "smooth") + 
  # Add correlation coefficient
  stat_cor(
    #method = "pearson" 
    method = "spearman"
    #,label.x = 3, label.y = 30
              ) + 
  ggtitle(
    "Fronto-parietal Mask"
    #"Sensorimotor Mask"
    ) +
  theme(plot.title = element_text(hjust = 0.5))
  
