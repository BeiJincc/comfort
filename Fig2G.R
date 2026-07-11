############################################################
# Figure 2G
# Clone size distribution across response groups
############################################################

library(ggplot2)

#-----------------------------------------------------------
# Prepare proportion table
#-----------------------------------------------------------

proportion_df <- as.data.frame(
  prop.table(
    table(
      bcr_all$response_group,
      bcr_all$seq_count
    )
  )
)

colnames(proportion_df) <- c(
  "response_group",
  "clone_size",
  "proportion"
)

group_colors <- colorRampPalette(
  c(
    "#619FBDFF",
    "#8CB9CFFF",
    "#B7D3E1FF",
    "#DEE9ECFF",
    "#faeebf",
    "#FFCE4EFF",
    "#FBA23CFF",
    "#FF7F00",
    "#E41A1C"
  )
)(10)

p_clone_proportion <- ggplot(
  proportion_df,
  aes(
    x = response_group,
    y = proportion,
    fill = clone_size
  )
) +
  geom_bar(stat="identity",position="fill") +
  theme_bw() +
  scale_fill_manual(values = group_colors) +
  guides(fill=guide_legend(title=NULL)) +
  labs(
    x=NULL,
    y="Proportion",
    title="Clone Size Distribution"
  ) +
  theme(
    axis.text=element_text(color="black"),
    axis.title=element_text(color="black"),
    axis.text.x=element_text(angle=45,hjust=1),
    panel.grid=element_blank(),
    panel.border=element_rect(color="black"),
    axis.line=element_line(color="black")
  )

ggsave(
  "Figure2G_clone_proportion.pdf",
  p_clone_proportion,
  width=3,
  height=5
)
