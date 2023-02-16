library(tidyr)
library(lmerTest)
library(lattice)

### Linear Mixed Effects Model of the R-baclofen data

## Load and arrange the data
d <- read.csv("ASC_TD_Pla_15mg_30mg_wDC_Z_hyper.csv", header = TRUE, row.names = 1)
d$id <- 1:nrow(d)
colnames(d) <- c("Group","Baseline", "x15mg", "x30mg")
d2 <- gather(data = d, key = "Condition", value = "wDC", Baseline, x15mg, x30mg)
colnames(d2) <- c("Group", "id", "Condition", "wDC")
d2[ d2 == "Baseline" ] <- 0
d2[ d2 == "x15mg" ] <- 15
d2[ d2 == "x30mg" ] <- 30
d2$Condition <- as.double(d2$Condition)

## Linear Mixed Effects Model
# https://ademos.people.uic.edu/Chapter18.html
# https://benwhalley.github.io/just-enough-r/extending-traditional-rm-anova.html
wDC.lmer <- lmer(wDC ~ Group:Condition + (1|id), data=d2, na.action=na.exclude)

# Check the summary results of the model
summary(wDC.lmer)
# anova(wDC.lmer)

# Check the random effects
ranef(wDC.lmer)

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

##