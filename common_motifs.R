common_motifs <- read.delim('./chipseq/commonmotifs.txt', na.strings = c("", NA))

common_motifs[c]
colnames(common_motifs)

sub_comm <- common_motifs[,c(22:71)]
sub_comm_matrix <- apply(sub_comm, 2, function(x) ifelse(is.na(x), 0, 1))

rownames(sub_comm_matrix) <- common_motifs[[1]]

colnames(sub_comm_matrix) <- lapply(
  colnames(sub_comm_matrix),
  function(x){unlist(strsplit(x, split="\\."))[2]})

colnames(sub_comm_matrix)

#HEATMAP FOR CLUSTERS

cluster_comm <- get_clusters(sub_comm_matrix)

heatmap(cluster_comm, Rowv = NA, Colv = NA)

df

round((df[1,]/25000 * 100/ df[2,]/13000 * 100))
