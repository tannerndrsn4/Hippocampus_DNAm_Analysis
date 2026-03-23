library(dplyr)

cluster_map <- c(
  "Cluster1" = "Cluster2",
  "Cluster2" = "Cluster4",
  "Cluster3" = "Cluster1",
  "Cluster4" = "Cluster3"
)

# Apply to all datasets
TE_candidates <- TE_candidates %>%
  mutate(DNAm_Cluster = cluster_map[DNAm_Cluster])

Enhancer_candidates <- Enhancer_candidates %>%
  mutate(DNAm_Cluster = cluster_map[DNAm_Cluster])

genic_candidates <- genic_candidates %>%
  mutate(DNAm_Cluster = cluster_map[DNAm_Cluster])

write.csv(TE_candidates,
          "TE_candidates_corrected_clusters.csv",
          row.names = FALSE)

write.csv(Enhancer_candidates,
          "Enhancer_candidates_corrected_clusters.csv",
          row.names = FALSE)

write.csv(genic_candidates,
          "Genic_candidates_corrected_clusters.csv",
          row.names = FALSE)

TE_df <- TE_candidates %>%
  select(RegionID, DNAm_Cluster) %>%
  mutate(Type = "TE")

Enhancer_df <- Enhancer_candidates %>%
  select(RegionID, DNAm_Cluster) %>%
  mutate(Type = "Enhancer")

intergenic_combined <- bind_rows(TE_df, Enhancer_df)

intergenic_summary <- intergenic_combined %>%
  group_by(RegionID, DNAm_Cluster) %>%
  summarise(
    Type = case_when(
      all(Type == "TE") ~ "TE",
      all(Type == "Enhancer") ~ "Enhancer",
      TRUE ~ "Both"
    ),
    .groups = "drop"
  )

genic_summary <- genic_candidates %>%
  group_by(feature) %>%
  summarise(n = dplyr::n()) %>%
  mutate(Group = "Genic",
         Category = feature)

intergenic_feature_summary <- intergenic_summary %>%
  group_by(Type) %>%
  summarise(n = dplyr::n()) %>%
  mutate(Group = "Intergenic",
         Category = Type)

feature_plot_df <- bind_rows(genic_summary, intergenic_feature_summary)

library(ggplot2)

feature_levels <- c(
  "Promoter",
  "Intron",
  "Exon",
  "UTR",
  "Enhancer",
  "TE"
)

feature_plot_df$Category <- factor(feature_plot_df$Category,
                                   levels = feature_levels)

feature_colors <- c(
  "Promoter" = "#0072B2",
  "Intron"   = "#009E73",
  "Exon"     = "#E69F00",
  "UTR"      = "#CC79A7",
  "Enhancer" = "#D55E00",
  "TE"       = "#56B4E9"
)


plot1 <- ggplot(feature_plot_df,
                aes(x = Group, y = n, fill = Category)) +
  geom_bar(stat = "identity",
           width = 0.55,
           color = "black",
           linewidth = 0.4) +   # outline thickness
  scale_fill_manual(values = feature_colors) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "Number of candidate regions",
    fill = "Feature"
  )

plot1

ggsave("Candidate_feature_breakdown.pdf",
       plot1,
       width = 6,
       height = 4)


genic_cluster_summary <- genic_candidates %>%
  group_by(DNAm_Cluster) %>%
  summarise(n = dplyr::n()) %>%
  mutate(Group = "Genic")

intergenic_cluster_summary <- intergenic_summary %>%
  group_by(DNAm_Cluster) %>%
  summarise(n = dplyr::n()) %>%
  mutate(Group = "Intergenic")

cluster_plot_df <- bind_rows(genic_cluster_summary,
                             intergenic_cluster_summary)


cluster_colors <- c(
  "Cluster1" = "#332288",  # deep indigo
  "Cluster2" = "#AA4499",  # plum
  "Cluster3" = "#882255",  # wine
  "Cluster4" = "#DDCC77"   # olive
)

plot2 <- ggplot(cluster_plot_df,
                aes(x = Group, y = n, fill = DNAm_Cluster)) +
  geom_bar(stat = "identity",
           width = 0.55,
           color = "black",
           linewidth = 0.4) +
  scale_fill_manual(values = cluster_colors) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "Number of candidate regions",
    fill = "DNAm cluster"
  )

plot2

ggsave("Candidate_cluster_breakdown.pdf",
       plot2,
       width = 6,
       height = 4)

library(tidyverse)

TE_long <- TE_candidates %>%
  separate_rows(TE_targets, sep = ";\\s*") %>%
  filter(!is.na(TE_targets))

TE_counts <- TE_long %>%
  count(TE_targets, sort = TRUE) %>%
  slice_head(n = 20)   # top 20 TEs

plot3 <- ggplot(TE_counts,
       aes(x = reorder(TE_targets, n), y = n)) +
  geom_col(fill = "#56B4E9", color = "black") +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    x = "Transposable element",
    y = "Number of intersecting DNAm regions"
  )
plot3

ggsave("TEsIntersectingDNAmRegions.pdf",
       plot3,
       width = 6,
       height = 4)

TE_cluster_summary <- TE_candidates %>%
  group_by(DNAm_Cluster) %>%
  summarise(n = dplyr::n()) %>%
  mutate(Group = "TE")

Enhancer_cluster_summary <- Enhancer_candidates %>%
  group_by(DNAm_Cluster) %>%
  summarise(n = dplyr::n()) %>%
  mutate(Group = "Enhancer")

cluster_TE_enhancer_df <- bind_rows(
  TE_cluster_summary,
  Enhancer_cluster_summary
)

plot4 <- ggplot(cluster_TE_enhancer_df,
       aes(x = Group, y = n, fill = DNAm_Cluster)) +
  geom_bar(stat = "identity",
           width = 0.55,
           color = "black") +
  scale_fill_manual(values = cluster_colors) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "Number of intergenic candidate regions",
    fill = "DNAm cluster"
  )

plot4

ggsave("TE_Enhancer_DNAmCluster_Breakdown.pdf",
       plot4,
       width = 6,
       height = 4)



library(ggVennDiagram)
library(rlang)

# ---- create dummy elements ----
genic_unique <- paste0("G", 1:467)
enhancer_unique <- paste0("E", 1:188)
shared <- paste0("S", 1:25)

gene_list <- list(
  "Genic targets"    = c(genic_unique, shared),
  "Enhancer targets" = c(enhancer_unique, shared)
)

p <- ggVennDiagram(gene_list, label_alpha = 0) +
  scale_fill_gradient(low = "#FFFFFF", high = "#D55E00") +
  theme(
    legend.position = "none",
    text = element_text(size = 14)
  )

# ---- REMOVE the set-name text layer (the one causing the cutoff) ----
for (i in seq_along(p$layers)) {
  lyr <- p$layers[[i]]
  
  # only touch layers that draw labels
  if (!is.null(lyr$mapping$label)) {
    lab_expr <- rlang::as_label(lyr$mapping$label)
    
    # ggVennDiagram uses a label aesthetic for set/category names (usually "name")
    if (grepl("\\bname\\b", lab_expr)) {
      p$layers[[i]]$aes_params$alpha <- 0  # hide those labels
    }
  }
}

p

ggsave("Genic_EnhancerOverlap_Venn.pdf",
       p,
       width = 6,
       height = 4)


# ---- create dummy elements ----
genic_unique <- paste0("G", 1:2592)
enhancer_unique <- paste0("E", 1:126)
shared <- paste0("S", 1:87)

gene_list <- list(
  "Genic targets"    = c(genic_unique, shared),
  "Enhancer targets" = c(enhancer_unique, shared)
)

p <- ggVennDiagram(gene_list, label_alpha = 0) +
  scale_fill_gradient(low = "#3abeab", high = "#D55E00") +
  theme(
    legend.position = "none",
    text = element_text(size = 14)
  )

# ---- REMOVE the set-name text layer (the one causing the cutoff) ----
for (i in seq_along(p$layers)) {
  lyr <- p$layers[[i]]
  
  # only touch layers that draw labels
  if (!is.null(lyr$mapping$label)) {
    lab_expr <- rlang::as_label(lyr$mapping$label)
    
    # ggVennDiagram uses a label aesthetic for set/category names (usually "name")
    if (grepl("\\bname\\b", lab_expr)) {
      p$layers[[i]]$aes_params$alpha <- 0  # hide those labels
    }
  }
}

p

ggsave("DEGs_EnhancerOverlap_Venn.pdf",
       p,
       width = 6,
       height = 4)