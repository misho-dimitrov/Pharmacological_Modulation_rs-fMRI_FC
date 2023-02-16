library(tidyr)
library(ggplot2)
library(ggpubr)

## Load and arrange the data
d <- read.csv("ASC_Pla_15mg_30mg_wDC_Z_hyper.csv", header = TRUE, row.names = 1)
n_asc <- nrow(d)
asc_subject_n <- paste("Subject (n=", n_asc, ")", sep="")
d$id <- rownames(d)
d2 <- gather(data = d, key = "Condition", value = "wDC", Placebo, X15mg, X30mg)

d3 <- read.csv("TD_Pla_15mg_30mg_wDC_Z_hyper.csv", header = TRUE, row.names = 1)
n_td <- nrow(d3)
td_subject_n <- paste("Subject (n=", n_td, ")", sep="")
d3$id <- rownames(d3)
d4 <- gather(data = d3, key = "Condition", value = "wDC", Placebo, X15mg, X30mg)

# Set min and max values
asc_min <- min(d2$wDC, na.rm = TRUE)
td_min <- min(d4$wDC, na.rm = TRUE)
overall_min <- pmin(asc_min, td_min)

asc_max <- max(d2$wDC, na.rm = TRUE)
td_max <- max(d4$wDC, na.rm = TRUE)
overall_max <- pmax(asc_max, td_max)

## Plot
ASC <- ggplot(d2, aes(x = Condition, y = wDC, group = id, colour = as.factor(id))) + 
  geom_point(na.rm = TRUE)+ 
  geom_line(na.rm = TRUE)+
  stat_summary(aes(group = 1), fun = median, geom = "point", shape = 18, size = 6, na.rm = TRUE)+
  stat_summary(aes(group = 1), fun = median, geom = "line",size = 1, na.rm = TRUE) +
  labs(title = "ASC", color = asc_subject_n) + ylim(overall_min, overall_max)

TD <- ggplot(d4, aes(x = Condition, y = wDC, group = id, colour = as.factor(id))) + 
  geom_point(na.rm = TRUE)+ 
  geom_line(na.rm = TRUE)+
  stat_summary(aes(group = 1), fun = median, geom = "point", shape = 18, size = 6, na.rm = TRUE)+
  stat_summary(aes(group = 1), fun = median, geom = "line",size = 1, na.rm = TRUE) +
  labs(title = "TD", color = td_subject_n) + ylim(overall_min, overall_max) 

final_figure <- ggarrange(ASC, TD, ncol = 2, nrow = 1)

annotate_figure(final_figure, top = text_grob("Differences in Mean Fronto-parietal Weighted Degree Centrality", 
                                      #color = "red", face = "bold", size = 14
                                      ))
