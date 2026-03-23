library(ggplot2)
library(dplyr)

var_results <- read.csv(file = "Desktop/PhD brain evolution/DNA methylation/Variance Analysis/topVar_Age_Sex_SV1-3_FDR0.1_results.csv")

# Add useful columns
var_results <- var_results %>%
  mutate(
    Direction = ifelse(LogVarRatio > 0, "Increased", "Decreased"),
    negLogFDR = -log10(Age_FDR)
  )

# Volcano plot with vertical line at x = 0
volcano_plot <- ggplot(var_results, aes(x = LogVarRatio, y = negLogFDR, color = Direction)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("Increased" = "#fd8d3c", "Decreased" = "#756bb1"),
    labels = c("Increased" = "Increased across age",
               "Decreased" = "Decreased across age"),
    guide = guide_legend(override.aes = list(size = 4))  # Larger legend dots
  ) +
  labs(
    x = "LogVarianceRatio",
    y = "-log10(FDR)",
    color = "Variance Change"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

volcano_plot

# Save as PDF
ggsave(
  filename = "volcano_variance_plot.pdf",
  plot = volcano_plot,
  path = "~/Desktop/PhD brain evolution/DNA methylation/Variance Analysis",
  width = 7, height = 5, units = "in"
)
