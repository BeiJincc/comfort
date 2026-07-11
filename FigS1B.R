############################################################
# Figure S1B
# GO Biological Process enrichment
############################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(forcats)
library(readxl)
library(Seurat)
library(Matrix)

table(seu_plasma$response_group)
# Identify differentially expressed genes
deg_results <- Finddeg_resultss(seu_plasma,group.by = "response_group",ident.1 ="Non_responders" ,logfc.threshold = 0.25)
write.csv(deg_results,"Plasma_Non_responders_DEG.csv")

rownames(deg_results)
data_upregulated<-rownames(deg_results)[deg_results$avg_log2FC > 0]
data_downregulated<-rownames(deg_results)[deg_results$avg_log2FC < 0]

# Convert gene list to character vectors
data_upregulated<-as.character(data_upregulated)
data_downregulated<-as.character(data_downregulated)

# Convert gene symbols to Entrez IDs
deg_up_entrez<- bitr(data_upregulated,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
head(deg_up_entrez,2)

deg_down_entrez<- bitr(data_downregulated,fromType = "SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")
head(deg_down_entrez,2)

# groupGO enrichment analysis
#ggo<-groupGO(gene=test1$ENTREZID,OrgDb = org.Hs.eg.db,ont = "CC",level=3,readable = TRUE)

# GO enrichment analysis
# Comment
#ego_ALL<-enrichGO(gene=test1$ENTREZID,
#OrgDb = org.Hs.eg.db,
# Comment
# Comment
# Comment
#qvalueCutoff = 0.05,
# Comment
#head(ego_ALL,2)

# Biological Process enrichment
go_bp_up<-enrichGO(gene=deg_up_entrez$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE) 
head(go_bp_up,2)
go_bp_down<-enrichGO(gene=deg_down_entrez$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE) 
head(go_bp_down,2)

write.csv(summary(go_bp_up),"Non_responders_up_enrich_0.05.csv",row.names = FALSE)
write.csv(summary(go_bp_down),"Non_responders_down_enrich_0.05.csv",row.names = FALSE)

# Save figure
library(ggplot2)
# Dot plot
p_go_bp <- dotplot(go_bp_up,showCategory = 15,title="EnrichmentGO_upregulated_BP_dot")
print(p_go_bp)
p_go_bp$data
ggsave("Non_responders_up_enrich_0.05.pdf",width = 8,height = 10)