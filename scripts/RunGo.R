#
library(clusterProfiler)
library(org.Mm.eg.db)

#read DEG's from output folder

ego_kw <- enrichGO(gene          = rownames(KW),
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                keyType = "SYMBOL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

barplot(ego_kw, showCategory = 20)
head(ego_kw)
