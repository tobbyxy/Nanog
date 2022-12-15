#script to get differentially expressed genes.

#read in meta data

targets <- read.csv("./SraRunTable (1).txt")

path <- file.path(targets$sample, "./abundance.tsv")
all(file.exists(path))

#read in transcript data 
Tx <- transcripts(EnsDb.Mmusculus.v79, columns=c("tx_id", "gene_name"))

Tx <- as.tibble(Tx)

#format column name
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")

#use Txi import to read count data from kallisto
Txi_gene <- tximport::tximport(path, type="kallisto",
                               tx2gene = Tx,
                               txOut = FALSE,
                               countsFromAbundance = "lengthScaledTPM",
                               ignoreTxVersion = T)

Txi_gene_no <- tximport::tximport(path, type="kallisto",
                                  tx2gene = Tx,
                                  txOut = FALSE,
                                  countsFromAbundance = "no",
                                  ignoreTxVersion = T)

#transcript level

myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts
#colnames(myCounts) < sampleLabels


sampleLabels <- targets$sample

myDGEList <- DGEList(myCounts)

#save Deg list 
save(myDGEList, file="myDGEList")


#get counts in millions
cpm <- cpm(myDGEList)


log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames= "geneID")

colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <-pivot_longer(log2.cpm.df, cols=NanogKO_0h:WT_24h, names_to = "samples", values_to = "expression")


#filter Degs for only 
keepers <- rowSums(cpm>1)>=3

myDGEList.filtered <- myDGEList[keepers,]
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)

log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames= "geneID")

colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <-pivot_longer(log2.cpm.filtered.df, cols=NanogKO_0h:WT_24h, names_to = "samples", values_to = "expression")


#normalize filtered data

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")

log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)

log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames= "geneID")

colnames(log2.cpm.filtered.norm.df) <- c("geneID", target$id)
colnames(log2.cpm.filtered.norm) <- c("geneID", target$id)
dimnames(log2.cpm.filtered.norm)[[2]] <- target$id



#get design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

#estimate mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot=TRUE)
fit <- lmFit(v.DEGList.filtered.norm,design)

#comparison between timepoints
contrast.matrix <- makeContrasts(k240= NanogKO_24h - NanogKO_0h,
                                 k120= NanogKO_12h - NanogKO_0h,
                                 w240=WT_24h -WT_0h,
                                 w120=WT_12h - WT_0h,
                                 levels=design)
contrast.matrix
fits <- contrasts.fit(fit, contrast.matrix)
ebfit <- eBayes(fits)
ebfit

mytophits_K240 <- topTable(ebfit, adjust="BH", coef = 1, number =40000, sort.by="logFC")
mytophits_W240 <- topTable(ebfit, adjust="BH", coef = 3, number =40000, sort.by="logFC")
mytophits_W120 <- topTable(ebfit, adjust="BH", coef = 4, number = 40000, sort.by = "logFC")
mytophits_K120 <- topTable(ebfit, adjust="BH", coef = 2, number = 40000, sort.by = "logFC")

#identify Degs
results <- decideTests(ebfit, method = "global", adjust.method = "BH", p.value = 0.05, lfc=1)
summary(results)
vennDiagram(results[,c(1,3)],include = "up")

#get Degs that meets conditions
#if gene is upregulated in knockout and don't change expression in wid type - downregulated
#if gene is upregulated in wild type and don't change expression in knockout - upregulated
downreg_genes <- v.DEGList.filtered.norm$E[results[,1] == 1 & results[,3] == 0,]
upreg_genes <- v.DEGList.filtered.norm$E[results[,1] == 0 & results[,3] == 1,]

all_degs <- v.DEGList.filtered.norm$E[results[,1] == 1 | results[,1] == -1 | results[,2] == 1 | results[,2] == -1,]

head(all_degs)

colnames(diffgenes) <- targets$id
colnames(upgenes) <- targets$id

dim(upgenes)

#save list for downstream analysis
write.table(rownames(downreg_gene), "downgenes.txt", quote = F, col.names = F, row.names = F)
write.table(rownames(upreg_genes), "upgenes.txt", quote =F, col.names = F, row.names = F)


