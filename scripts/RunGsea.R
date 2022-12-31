#
library(msigdbr)

msigdbr_show_species()

ms_gsea <- msigdbr(species = "Mus musculus")
ms_gsea %>% distinct(gs_cat, gs_subcat) %>% 
  arrange(gs_cat, gs_subcat)

ms_gsea_c2 <- msigdbr(species = "Mus musculus", category = "C2") %>% 
  dplyr::select(gs_name, gene_symbol)


#construct a named vector

kw.gsea <- KW$log2FoldChange
names(kw.gsea) <- as.character(rownames(KW))
kw.gsea <- sort(kw.gsea, decreasing = T)

#run GSEA using the GSEA function from cluster profiler

kw.gsea.res <- GSEA(kw.gsea, TERM2GENE = ms_gsea_c2, verbose = F)
kw.gsea.df <- kw.gsea.res@result
datata


kw.ego <- gseGO(geneList     = kw.gsea,
                keyType = "SYMBOL",
              OrgDb        = org.Mm.eg.db,
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
dotplot(kw.ego)
head(kw.gsea.df[order(kw.gsea.df$enrichmentScore, decreasing = T),])

