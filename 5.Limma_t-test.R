# modify the column numbers in lines 13 and 14 to reflect the specific columns you wish to compare (count)
# In line 31, update the number of replicates to match the number of input columns specified in lines 13 and 14.

# Load the necessary libraries
library(limma)
library(ggplot2)
library(ggrepel)

# Read the dataset
data <- read.csv("file_for_R_script.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Extract quantitative values for each condition
Light_values <- data[, c(7,9,11)]  # Load columns with light instensity
Heavy_values <- data[, c(8,10,12)]  # Load columns with heavy instensity

# Log2 transform the data
Light_values <- log2(Light_values)  
Heavy_values <- log2(Heavy_values)  

# Set the desired column names for each group
colnames(Light_values) <- paste("Light", seq_len(ncol(Light_values)), sep = "_")
colnames(Heavy_values) <- paste("Heavy", seq_len(ncol(Heavy_values)), sep = "_")

# Combine all data into a single matrix with updated column names
all_data <- cbind(Light_values, Heavy_values)

# Convert relevant columns in all_data to numeric
all_data_numeric <- apply(all_data, 2, as.numeric)

# Define the number of replicates for each condition
num_replicates <- c(3, 3)  # For CTL, STR respectively

# Create factor levels based on the number of replicates
factor_levels <- rep(c("Light", "Heavy"), times = num_replicates)

# Define the design matrix
design <- model.matrix(~0 + factor(factor_levels))

# Set appropriate column names for the design matrix
colnames(design) <- c("Light","Heavy")

# Fit the linear model
fit <- lmFit(all_data_numeric, design)
head(coef(fit))

# Define contrasts
contrasts <- makeContrasts(
  Light_vs_Heavy = Heavy - Light,
  levels = colnames(coef(fit))  # Use the moderation object for levels
)
contrasts

# Perform contrasts
contrast_fit <- contrasts.fit(fit, contrasts)

# Perform empirical Bayes moderation
moderation <- eBayes(contrast_fit)

# Extract moderated t-statistics and p-values for CTL_vs_STR
moderated_t_Light_vs_Heavy <- moderation$coefficients[, "Light_vs_Heavy"]
p_values_Light_vs_Heavy <- moderation$p.value[,"Light_vs_Heavy"]
# choose how and if you want to adjust your p value
adj_p_values_Light_vs_Heavy <- p.adjust(p_values_Light_vs_Heavy, method = "fdr")
#adj_p_values_WT_vs_mutant <- p.adjust(p_values_WT_vs_mutant, method = "bonferroni")
#adj_p_values_WT_vs_mutant <- p.adjust(p_values_WT_vs_mutant, method = "holm")

# Add the extracted statistics to your data matrix
data$moderated_t_Light_vs_Heavy <- moderated_t_Light_vs_Heavy
data$p_values_Light_vs_Heavy <- p_values_Light_vs_Heavy
data$adj_p_values_Light_vs_Heavy <- adj_p_values_Light_vs_Heavy

# Convert CTL_values and STR_values to numeric format
Light_values_numeric <- apply(Light_values, 2, as.numeric)
Heavy_values_numeric <- apply(Heavy_values, 2, as.numeric)

# Check for non-numeric values in CTL_values and STR_values
if (any(is.na(Light_values_numeric)) || any(is.na(Heavy_values_numeric))) {
  print("Warning: There are non-numeric values in your data.")
  # You might need to handle missing or non-numeric values appropriately.
}

# Create a dataframe with average values
average_data <- data.frame(
  Average_Light = rowMeans(Light_values_numeric, na.rm = TRUE),
  Average_Heavy = rowMeans(Heavy_values_numeric, na.rm = TRUE))

# Append average_data to data
data_final <- cbind(data, average_data)

# Calculate the logFC for each comparison
data_final$Heavy_over_Light_logFC <- data_final$Average_Heavy - data_final$Average_Light

write.table(data_final, file = "Output_R.tsv", sep = "\t", row.names = FALSE)

