############################################################
# Figure S1D
# IGHV usage in clonally expanded plasma cells
############################################################

# Separate clonally expanded and unique plasma cells
bcr_plasma_clone<-subset(bcr_plasma,clone_group %in% "Clone")
bcr_plasma_unique<-subset(bcr_plasma,clone_group %in% "Unique")

#########response_group clone IGHV#####
library(dplyr)
library(tidyr)
library(ggplot2)

# Build contingency table
summary_df <- as.data.frame.matrix(table(bcr_plasma_clone$v_call_10x, bcr_plasma_clone$response_group))
summary_df$v_call_10x <- rownames(summary_df)

# Calculate proportions and log2 fold change
summary_df <- summary_df %>%
  mutate(
    Responders_prop = Responders / sum(Responders),
    Non_responders_prop = Non_responders / sum(Non_responders),
    log2FC = log2(
      (Non_responders_prop + 1e-6) /
      (Responders_prop + 1e-6)
    ),
    total_count = Responders + Non_responders
  )

# Statistical test
calculate_p_value <- function(x) {

  mat <- matrix(
    c(
      as.numeric(x["Responders"]),
      as.numeric(x["Non_responders"]),
      sum(summary_df$Responders) - as.numeric(x["Responders"]),
      sum(summary_df$Non_responders) - as.numeric(x["Non_responders"])
    ),
    nrow = 2
  )

  if (any(chisq.test(mat)$expected < 5)) {
    fisher.test(mat)$p.value
  } else {
    chisq.test(mat)$p.value
  }
}
summary_df$p_value <- apply(summary_df, 1, calculate_p_value)

# Add significance labels
summary_df$signif_label <- cut(summary_df$p_value,
                        breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                        labels = c("***", "**", "*", "ns"))

# Convert to long format
plot_data <- summary_df %>%
  select(
    v_call_10x,
    Responders_prop,
    Non_responders_prop,
    log2FC,
    p_value,
    signif_label,
    total_count
  ) %>%
  pivot_longer(
    cols = c(
      Responders_prop,
      Non_responders_prop
    ),
    names_to = "group",
    values_to = "proportion"
  ) %>%
  mutate(
    group = ifelse(
      group == "Responders_prop",
      "Responders",
      "Non_responders"
    ),
    direction = ifelse(
      group == "Responders",
      "up",
      "down"
    ),
    proportion_plot = ifelse(
      direction == "up",
      proportion,
      -proportion
    )
  )

# Order genes by total count
plot_data$v_call_10x <- factor(plot_data$v_call_10x,
                             levels = summary_df %>% arrange(desc(total_count)) %>% pull(v_call_10x))
print(plot_data, n = Inf) 
# Generate plot
ggplot(plot_data, aes(x = v_call_10x, y = proportion_plot, fill = group)) +
  geom_bar(stat = "identity", position = "identity", width = 0.8) +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(
  values = c(
    "Responders" = "#619FBDFF",
    "Non_responders" = "#E41A1C"
  )
) +
  labs(x = "v_call_10x", y = "Proportion", fill = "Group") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.grid = element_blank()
  ) +
  # Add significance labels
  geom_text(data = summary_df,
            aes(x = v_call_10x, y = 0.23, label = signif_label),
            inherit.aes = FALSE, size = 4, vjust = 0, color = "black") +
  # Add log2 fold-change labels
  geom_text(data = summary_df,
            aes(x = v_call_10x, y = 0.2, label = paste0("", round(log2FC, 2))),
            inherit.aes = FALSE, size = 3, vjust = 0, color = "black")

ggsave("pc_clone_IGHV.pdf", width = 18, height = 5)