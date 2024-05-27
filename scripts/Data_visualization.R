counts_data <- read.table('../../Downloads/counts_data.txt')
head(counts_data)


# read in sample info
colData <- read.table('../../Downloads/sample_info.txt')
colData
colnames(colData) <- "group"
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ group)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

counts(dds)
vsd <- vst(dds, blind=FALSE)
vsd
head(assay((vsd)))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colData$group
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("./Distance_matrix.pdf", width = 8, height= 8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
