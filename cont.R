
#script to get heat map of motifs

#read in motifs files
#merged means motifs found from both time points

all_motifs_00 <- read.delim('./Motif_matrix/all_motif_00.txt', na.strings = c("", "NA"))
all_motifs_24 <- read.delim('./Motif_matrix/all_motif_24.txt', na.strings = c("", "NA"))
merged_mtf_00 <- read.delim('./Motif_matrix/0hrmerged_motif.txt', na.strings = c("", "NA"))
merged_mtf_24 <- read.delim('./Motif_matrix/24hrmerged_motif.txt', na.strings = c("", "NA"))
common_motifs <- read.delim('./Motif_matrix/commonmotifs.txt', na.strings = c("", NA))


#get motif columns
#for all motif column 22-46



#for merged motifs column 22-71

all_motif_list <- list(all_motifs_00, all_motifs_24)
merged_motif_list <- list(merged_mtf_00, merged_mtf_24, common_motifs)

all_mtf <- lapply(all_motif_list, function(x) x[, c(22:46)])
merged_mtf <- lapply(merged_motif_list, function(x) x[, c(22:71)])


#merge files into a list
all_list <- list(all_motifs_00, all_motifs_24, common_motifs)


#get peak id of all file

peak_ids <- lapply(all_list, function (x) x[[1]])
peak_ids[[2]]

#get matrix files
all_matrix_list <- c(all_mtf,merged_mtf)

#using motif distance, if motif occurs 1 else 0

all_motif_matrix <- lapply(all_matrix_list, function(x) ifelse(is.na(x), 0, 1))

#change column names 

new_col_name <- lapply(all_motif_matrix, function (x){
  matrix_col_name(x)
})

for (i in seq_along(all_motif_matrix)){
  colnames(all_motif_matrix[[i]]) <- new_col_name[[i]]
}

#change rownames of files
rownames(all_motif_matrix[[1]]) <- peak_ids[[1]]
rownames(all_motif_matrix[[2]]) <- peak_ids[[2]]
rownames(all_motif_matrix[[3]]) <- peak_ids[[1]]
rownames(all_motif_matrix[[4]]) <- peak_ids[[2]]
rownames(all_motif_matrix[[5]]) <- peak_ids[[3]]


#get clusters among motifs
get_clusters(all_motif_matrix)
all_motif_matrix[c(3:5)]

cluster_list <- lapply(all_motif_matrix[c(3:5)], function(x) get_clusters(x))

#get heatmaps
heatmap(cluster_list[[1]], Colv = NA, Rowv = NA )

lapply(cluster_list, function(x) pheatmap(x, annotation_col = sample_col,
                                          cluster_rows = TRUE))

sample_col <- read.table('./table.txt', col.names = c("cluster", "timepoint"))
row.names(sample_col) <- sample_col[,1]
sample_col[,1] <- NULL

tmp <- get_clusters(all_motif_matrix[[5]])
pheatmap(tmp, annotation_col = sample_col)

#7/27

#extract all numbers from columns

dm_00 <- apply(m_mtf_00, 2, function(x) str_extract(x,"^-?\\d+"))

#convert all columns to numeric
sapply(dm_00, class)

dm_00_df <- as.data.frame(dm_00)
str(dm_00_df)

dm_00_df[] <- lapply(dm_00_df, function(x) as.numeric(as.character(x)))
dm_00_df %>% summarise_if(is.numeric, max)

dm_00_df[is.na(dm_00_df)] <- 0
max(dm_00_df[,2])


tmp <- dm_00_df[[3]]

sapply(dm_00_df, function(x) rowSums(tmp, x))
rowSums(tmp)
apply(dm_00_df[, c(1:2), 1, sum])
dm_00_df[[3]] - tmp
