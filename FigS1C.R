############################################################
# Figure S1C
# Isotype distribution
############################################################

# Compare immunoglobulin isotypes
table(response_group_db$c_call)
proportion_df <- as.data.frame(prop.table(table(response_group_db$celltype_fine, response_group_db$c_call, response_group_db$response_group), margin = 3))

# Set column names
colnames(proportion_df) <- c("celltype_fine", "c_call", "response_group", "proportion")
group_colors <- brewer.pal(n = 10, name = "Paired")
#group_colors <- colorRampPalette(c( "#619FBDFF", "#8CB9CFFF", "#B7D3E1FF", "#DEE9ECFF",'#faeebf',"#FFCE4EFF","#FBA23CFF","#FF7F00", "#E41A1C"))(54)
# Plot stacked bar chart
ggplot(proportion_df, aes(x = response_group, y = proportion, fill = c_call)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~celltype_fine, ncol = 10) + 
  ggtitle("Isotype") +
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.2, 'cm'),  
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
  scale_fill_manual(values = group_colors)
ggsave("BCR_isotype_celltype_fine_response_group.pdf",width = 11,height = 3.5)