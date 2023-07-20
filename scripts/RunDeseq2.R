## script to perform differential gene expression analysis using DESeq2 package
# setwd("/home/tobi/projects/nanog/Nanog/data")

# load libraries
library(DESeq2)
library(tidyverse)


# Step 1: preparing count data ----------------

# read in counts data
counts_data <- read.table('./data/counts_data.txt')
head(counts_data)


# read in sample info
colData <- read.table('./data/sample_info.txt')
colData
colnames(colData) <- "group"

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ group)
dds
counts(dds)
group
# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

counts(dds)

#create design matrix for pairwise comparison

group <- as.factor(colData$group)
design <- model.matrix(~0+group)

colnames(design) <- levels(group)
design

contrast.matrix <- makeContrasts(k240= NanogKO_24h - NanogKO_0h,
                                 k120= NanogKO_12h - NanogKO_0h,
                                 w240=WT_24h -WT_0h,
                                 w120=WT_12h - WT_0h,
                                 kVsw=(NanogKO_24h - NanogKO_0h)-(WT_24h -WT_0h),
                                 levels=design)
contrast.matrix

# step 3 : Time series for Deseq
dds.tmp <- DESeqDataSetFromMatrix(countData=counts$counts,
                              colData=colData(ddsTxi),             # or colData=group,
                              design= ~condition + time + condition:time)
dds.tmp <- DESeq(dds.tmp, test="LRT", reduced = ~ condition + time)
resTC <- results(dds.tmp, contrast = c("condition", "KO", "WT"))

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
resKO_24vs0 <- results(dds,contrast = c("group", "NanogKO_24h", "NanogKO_0h"))

resKO_12vs0 <- results(dds , contrast = c("group", "NanogKO_12h", "NanogKO_0h"))                 
resWT_24vs0 <- results(dds, contrast = c("group", "WT_24h", "WT_0h"))
resWT_12vs0 <- results(dds, contrast = c("group", 'WT_12h', "WT_0h"))
kVsw <- results(dds, contrast= list(c("group_NanogKO_24h_vs_NanogKO_0h", "group_WT_24h_vs_NanogKO_0h")))
kVsw


# Explore Results ----------------
resultsNames(dds)
d

