# Step 1: Load the data and necessary libraries
predicted_ages <- readRDS(file = "Desktop/PhD brain evolution/DNA methylation/Rapamycin/predicted_ages_summary_rapa.rds")
library(ggplot2)

summary(lm(Chronological_Age ~ Group, data = predicted_ages))

# Step 2: Fit separate linear models for R and C groups to analyze predicted age vs. chronological age
model_R <- lm(lambda.min ~ Chronological_Age, data = predicted_ages, subset = Group == "R")
model_C <- lm(lambda.min ~ Chronological_Age, data = predicted_ages, subset = Group == "C")

# Step 3: Extract R^2 values for each group
R2_R <- summary(model_R)$adj.r.squared
R2_C <- summary(model_C)$adj.r.squared

# Step 4: Generate a scatter plot with regression lines, annotate with R^2 values and slopes
# Calculate slopes for each group
slope_R <- coef(lm(lambda.min ~ Chronological_Age, data = predicted_ages, subset = Group == "R"))["Chronological_Age"]
slope_C <- coef(lm(lambda.min ~ Chronological_Age, data = predicted_ages, subset = Group == "C"))["Chronological_Age"]

ggplot(predicted_ages, aes(x = Chronological_Age, y = lambda.min, color = Group)) +
  geom_point() +
  geom_smooth(method = "lm", aes(group = Group), se = FALSE) +
  labs(
    x = "Chronological Age",
    y = "Predicted Age",
    title = "Predicted Age vs. Chronological Age by Group"
  ) +
  annotate("text", 
           x = max(predicted_ages$Chronological_Age), 
           y = min(predicted_ages$lambda.min),
           label = paste0("R^2 (R): ", round(R2_R, 3), "\nSlope (R): ", round(slope_R, 3)),
           color = "blue", 
           hjust = 1, 
           size = 2.5) +
  annotate("text", 
           x = max(predicted_ages$Chronological_Age), 
           y = min(predicted_ages$lambda.min) + 1.5,
           label = paste0("R^2 (C): ", round(R2_C, 3), "\nSlope (C): ", round(slope_C, 3)),
           color = "red", 
           hjust = 1, 
           size = 2.5) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Step 5: Perform correlation tests for lambda.min vs. chronological age for each group
# Correlation test for Group R
cor_test_R <- cor.test(
  predicted_ages$lambda.min[predicted_ages$Group == "R"],
  predicted_ages$Chronological_Age[predicted_ages$Group == "R"]
)
cat("Correlation test for Group R:\n")
print(cor_test_R)

# Correlation test for Group C
cor_test_C <- cor.test(
  predicted_ages$lambda.min[predicted_ages$Group == "C"],
  predicted_ages$Chronological_Age[predicted_ages$Group == "C"]
)
cat("\nCorrelation test for Group C:\n")
print(cor_test_C)

# Step 6: Fit a linear model for all data to calculate residuals
residual_model <- lm(lambda.min ~ Chronological_Age, data = predicted_ages)

# Step 7: Add residuals from the model to the data frame
predicted_ages$residuals <- residuals(residual_model)

# Step 8: Perform a t-test to determine if residuals differ significantly by group
test_residuals <- t.test(residuals ~ Group, data = predicted_ages)

# Step 9: Summarize t-test results for residuals
cat("T-test results: Residual Epigenetic Age ~ Group\n")
print(test_residuals)

# Step 10: Fit a model to estimate the group effect on residual epigenetic age
group_model <- lm(residuals ~ Group, data = predicted_ages)
summary_group_model <- summary(group_model)

# Step 11: Summarize group effect results and extract beta value
cat("Linear regression results: Residual Epigenetic Age ~ Group\n")
print(summary_group_model)
beta_group <- coef(summary_group_model)["GroupR", "Estimate"]  # Change "GroupR" to the group of interest
cat("\nEpigenetic age acceleration in Group R compared to Group C:", beta_group, "years\n")

# Step 12: Visualize residual epigenetic age differences between groups
# Boxplot for residuals by group
ggplot(predicted_ages, aes(x = Group, y = residuals, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    x = "Treatment Group",
    y = "Residual Epigenetic Age (years)"
  ) +
  theme_minimal()

# Step 13: Alternatively, create a density plot for residuals by group
ggplot(predicted_ages, aes(x = residuals, fill = Group, color = Group)) +
  geom_density(alpha = 0.4) +
  labs(
    x = "Residual Epigenetic Age (years)",
    y = "Density"
  ) +
  theme_minimal()