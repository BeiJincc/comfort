############################################################
# Figure 2L
# IgG1-positive plasma cell clone expansion
############################################################

library(dplyr)
library(ggplot2)

#-----------------------------------------------------------
# Extract IgG1-positive plasma cells
#-----------------------------------------------------------
bcr_plasma_ighg1 <- subset(
  bcr_plasma,
  c_call %in% "IGHG1"
)

#-----------------------------------------------------------
# Define clone groups
#-----------------------------------------------------------
bcr_plasma_ighg1 <- bcr_plasma_ighg1 %>%
  mutate(
    clone_group = if_else(seq_count == 1, "Unique", "Clone")
  )

bcr_plasma_ighg1$clone_group <- factor(
  bcr_plasma_ighg1$clone_group,
  levels = c("Unique", "Clone")
)

#-----------------------------------------------------------
# Calculate clone proportions
#-----------------------------------------------------------
proportion_df <- as.data.frame(
  prop.table(
    table(
      bcr_plasma_ighg1$clone_group,
      bcr_plasma_ighg1$response_group
    ),
    margin = 2
  )
)

colnames(proportion_df) <- c(
  "clone_group",
  "response_group",
  "proportion"
)

#-----------------------------------------------------------
# Statistical test
#-----------------------------------------------------------
contingency_table <- table(
  bcr_plasma_ighg1$clone_group,
  bcr_plasma_ighg1$response_group
)

test_result <- if (all(chisq.test(contingency_table)$expected >= 5)) {
  chisq.test(contingency_table)
} else {
  fisher.test(contingency_table)
}

significance_label <- ifelse(
  test_result$p.value < 0.001, "***",
  ifelse(
    test_result$p.value < 0.01, "**",
    ifelse(
      test_result$p.value < 0.05, "*", "ns"
    )
  )
)

#-----------------------------------------------------------
# Plot
#-----------------------------------------------------------
clone_colors <- c(
  Unique = "#619FBDFF",
  Clone  = "#E41A1C"
)

p_clone_expansion <- ggplot(
  proportion_df,
  aes(
    x = response_group,
    y = proportion,
    fill = clone_group
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  geom_text(
    data = data.frame(
      x = 1.5,
      y = 1.05,
      label = significance_label
    ),
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5,
    color = "black"
  ) +
  scale_fill_manual(values = clone_colors) +
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black")
  ) +
  guides(fill = guide_legend(title = NULL)) +
  labs(
    x = NULL,
    y = "Proportion",
    title = "IgG1-positive Plasma Cell Clone Expansion"
  )

p_clone_expansion

#-----------------------------------------------------------
# Save figure
#-----------------------------------------------------------
ggsave(
  "Figure2L_IgG1_plasma_clone_expansion.pdf",
  plot = p_clone_expansion,
  width = 3,
  height = 3.5
)
