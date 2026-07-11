############################################################
# Figure 2H
# Clone expansion across B-cell subsets
############################################################

library(dplyr)
library(ggplot2)

bcr_all <- bcr_all %>%
  mutate(
    clone_group=if_else(seq_count==1,"Unique","Clone")
  )

bcr_all$clone_group <- factor(
  bcr_all$clone_group,
  levels=c("Unique","Clone")
)

proportion_df <- as.data.frame(
  prop.table(
    table(
      bcr_all$celltype_fine,
      bcr_all$clone_group,
      bcr_all$response_group
    ),
    margin=c(1,3)
  )
)

colnames(proportion_df) <- c(
  "celltype_fine",
  "clone_group",
  "response_group",
  "proportion"
)

test_results <- data.frame()

for(ct in unique(bcr_all$celltype_fine)){
  tbl <- table(
    bcr_all$clone_group[bcr_all$celltype_fine==ct],
    bcr_all$response_group[bcr_all$celltype_fine==ct]
  )
  if(all(dim(tbl)==c(2,2))){
    test <- if(all(chisq.test(tbl)$expected>=5)) chisq.test(tbl) else fisher.test(tbl)
    test_results <- rbind(
      test_results,
      data.frame(
        celltype_fine=ct,
        p_value=test$p.value,
        label=ifelse(test$p.value<0.001,"***",
          ifelse(test$p.value<0.01,"**",
          ifelse(test$p.value<0.05,"*","ns")))
      )
    )
  }
}

proportion_df <- left_join(proportion_df,test_results,by="celltype_fine")

clone_colors <- c(Unique="#619FBDFF",Clone="#E41A1C")

p_clone_expansion <- ggplot(
  proportion_df,
  aes(response_group, proportion, fill = clone_group)
) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~celltype_fine, ncol = 10) +
  geom_text(
    data = proportion_df %>%
      filter(
        clone_group == "Clone",
        response_group == "Responders"
        ) %>%
      distinct(celltype_fine, label),
    aes(
      x = 1.5,
      y = 1.05,
      label = label
    ),
    inherit.aes = FALSE,
    size = 5,
    color = "black"
  ) +
  scale_fill_manual(values = clone_colors) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black")
  ) +
  guides(fill = guide_legend(title = NULL)) +
  labs(x = NULL, y = "Proportion")


ggsave("Figure2H_clone_expansion.pdf",p_clone_expansion,width=11,height=3.5)
