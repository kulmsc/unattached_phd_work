args <- commandArgs(trailingOnly=TRUE)

gene_loc <- read.table("NCBI37.3.gene.loc", stringsAsFactors=F, header=F)
mag_res <- read.table("magma.genes.out", stringsAsFactors=F, header=T)

gene_loc <- gene_loc[gene_loc[,1] %in% mag_res[,1],]
gene_loc <- gene_loc[order(gene_loc[,1])[rank(mag_res[,1])],]
mag_res$gene <- gene_loc[,6]

write.table(mag_res, paste0("results/magma.", args[1], ".", args[2], ".", args[3], ".txt"), row.names = F, col.names = F, sep = "\t", quote = F)
