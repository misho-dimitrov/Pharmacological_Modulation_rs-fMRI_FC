library(foreign)

setwd('~/Downloads/PhD/Analysis/Clustering_Drug_Response/R-bac_clus/data/')

data <- read.spss("TRADa_Data_2.sav", to.data.frame=TRUE)

neurotypical_data <- data[grep("^0", data$ID), ]  # Subset data where ID starts with 0
autism_data <- data[grep("^1", data$ID), ]  # Subset data where ID starts with 1

#########################
#0-7: Normal range (no depression)
#8-16: Mild depression
#17-23: Moderate depression
#24-52: Severe depression

# Extract HRSDTotal values from the subset
HRSDTotal_NT_subset <- neurotypical_data$HRSDTotal
HRSDTotal_NT_subset <- as.numeric(as.character(HRSDTotal_NT_subset))

HRSDTotal_A_subset <- autism_data$HRSDTotal
HRSDTotal_A_subset <- as.numeric(as.character(HRSDTotal_A_subset))

# Plot the subsetted values against their index
#plot(HRSDTotal_subset, col = "blue", main = "Plot of HRSDTotal (Autism)")

# Histogram of HRSDTotal subset
#hist(HRSDTotal_subset, col = "lightblue", main = "Histogram of HRSDTotal (Neurotypical)")

# Boxplot of HRSDTotal subset
#boxplot(HRSDTotal_subset, col = "lightgreen", main = "Boxplot of HRSDTotal (Autism)")

# Calculate mean and standard deviation for group HRSDTotal_NT_subset, excluding NA values
mean_NT <- mean(HRSDTotal_NT_subset, na.rm = TRUE)
sd_NT <- sd(HRSDTotal_NT_subset, na.rm = TRUE)

# Calculate mean and standard deviation for group HRSDTotal_A_subset, excluding NA values
mean_A <- mean(HRSDTotal_A_subset, na.rm = TRUE)
sd_A <- sd(HRSDTotal_A_subset, na.rm = TRUE)

# Print mean and standard deviation for both groups
cat("Mean for HRSDTotal_NT_subset:", mean_NT, "\n")
cat("Standard Deviation for HRSDTotal_NT_subset:", sd_NT, "\n")

cat("\nMean for HRSDTotal_A_subset:", mean_A, "\n")
cat("Standard Deviation for HRSDTotal_A_subset:", sd_A, "\n")

# Run the Wilcoxon rank-sum test
wilcox.test(HRSDTotal_NT_subset, HRSDTotal_A_subset)

#########################
#0-17: Mild anxiety
#18-24: Moderate anxiety
#25-30: Severe anxiety
#31 or above: Very severe anxiety

# Extract HAMATotal values from the subset
HAMATotal_NT_subset <- neurotypical_data$HAMATotal
HAMATotal_NT_subset <- as.numeric(as.character(HAMATotal_NT_subset))

HAMATotal_A_subset <- autism_data$HAMATotal
HAMATotal_A_subset <- as.numeric(as.character(HAMATotal_A_subset))

# Plot the subsetted values against their index
#plot(HAMATotal_subset, col = "blue", main = "Plot of HAMATotal (Autism)")

# Histogram of HAMATotal subset
#hist(HAMATotal_subset, col = "lightblue", main = "Histogram of HAMATotal (Neurotypical)")

# Boxplot of HAMATotal subset
#boxplot(HAMATotal_subset, col = "lightgreen", main = "Boxplot of HAMATotal (Autism)")

# Calculate mean and standard deviation for group HAMATotal_NT_subset, excluding NA values
mean_HAMA_NT <- mean(HAMATotal_NT_subset, na.rm = TRUE)
sd_HAMA_NT <- sd(HAMATotal_NT_subset, na.rm = TRUE)

# Calculate mean and standard deviation for group HAMATotal_A_subset, excluding NA values
mean_HAMA_A <- mean(HAMATotal_A_subset, na.rm = TRUE)
sd_HAMA_A <- sd(HAMATotal_A_subset, na.rm = TRUE)

# Print mean and standard deviation for both groups
cat("Mean for HAMATotal_NT_subset:", mean_HAMA_NT, "\n")
cat("Standard Deviation for HAMATotal_NT_subset:", sd_HAMA_NT, "\n")

cat("\nMean for HAMATotal_A_subset:", mean_HAMA_A, "\n")
cat("Standard Deviation for HAMATotal_A_subset:", sd_HAMA_A, "\n")

# Run the Wilcoxon rank-sum test for HAMA
wilcox.test(HAMATotal_NT_subset, HAMATotal_A_subset)
