#script to get gene id 


library(biomaRt)



# retrieve gene info of our DEGs

ensembl = useMart ("ensembl", dataset="mmusculus_gene_ensembl")
retrievedGeneInfo_KW = getBM(attributes=c('external_gene_name','entrezgene_id', "hgnc_symbol"),
                          filters = c('external_gene_name'),
                          values = rownames(KW),
                          mart = ensembl, uniqueRows=TRUE)


