library(tidyr)
library(ggplot2)
library(ggpubr)

## Load and arrange the data
# Load data from autistic individuals
d <- read.csv("ASC_Pla_CBD_wDC_Z_hypo.csv", header = TRUE, row.names = 1)
n_asc <- nrow(d)
asc_subject_n <- paste("Subject (n=", n_asc, ")", sep="")
d$id <- rownames(d)
colnames(d) <- c("Baseline", "Drug", "id")

# Set trajectory colour rule 
asc_shift <- d$Baseline-d$Drug
asc_shift_colours <- ifelse(asc_shift > 0, 'red', 'blue')
asc_shift_colours_all <- c(asc_shift_colours, asc_shift_colours)

d2 <- gather(data = d, key = "Condition", value = "wDC", Baseline, Drug)

# Load data from neurotypical individuals
d3 <- read.csv("TD_Pla_CBD_wDC_Z_hypo.csv", header = TRUE, row.names = 1)
n_td <- nrow(d3)
td_subject_n <- paste("Subject (n=", n_td, ")", sep="")
d3$id <- rownames(d3)
colnames(d3) <- c("Baseline", "Drug", "id")

# Set trajectory colour rule 
td_shift <- d3$Baseline-d3$Drug
td_shift_colours <- ifelse(td_shift > 0, 'red', 'blue')
td_shift_colours_all <- c(td_shift_colours, td_shift_colours)

d4 <- gather(data = d3, key = "Condition", value = "wDC", Baseline, Drug)

# Set min and max values
asc_min <- min(d2$wDC, na.rm = TRUE)
td_min <- min(d4$wDC, na.rm = TRUE)
overall_min <- pmin(asc_min, td_min)

asc_max <- max(d2$wDC, na.rm = TRUE)
td_max <- max(d4$wDC, na.rm = TRUE)
overall_max <- pmax(asc_max, td_max)

## Plot
ASC <- ggplot(d2, aes(x = Condition, y = wDC, group = id, colour = as.factor(asc_shift_colours_all))) + 
  geom_point(na.rm = TRUE)+ 
  geom_line(na.rm = TRUE)+
  stat_summary(aes(group = 1), fun = median, geom = "point", shape = 18, size = 6, na.rm = TRUE)+
  stat_summary(aes(group = 1), fun = median, geom = "line",size = 1, na.rm = TRUE) +
  labs(title = "ASC", color = asc_subject_n) + 
  scale_color_manual(labels = c("increase", "decrease"), values = c("red", "blue")) +
  ylim(overall_min, overall_max)

TD <- ggplot(d4, aes(x = Condition, y = wDC, group = id, colour = as.factor(td_shift_colours_all))) + 
  geom_point(na.rm = TRUE)+ 
  geom_line(na.rm = TRUE)+
  stat_summary(aes(group = 1), fun = median, geom = "point", shape = 18, size = 6, na.rm = TRUE)+
  stat_summary(aes(group = 1), fun = median, geom = "line",size = 1, na.rm = TRUE) +
  labs(title = "TD", color = td_subject_n) + 
  scale_color_manual(labels = c("increase", "decrease"), values = c("red", "blue")) +
  ylim(overall_min, overall_max) 

final_figure <- ggarrange(TD, ASC, ncol = 2, nrow = 1)

annotate_figure(final_figure, top = text_grob("Differences in Mean Sensorimotor Weighted Degree Centrality", 
                                      #color = "red", face = "bold", size = 14
                                      ))
