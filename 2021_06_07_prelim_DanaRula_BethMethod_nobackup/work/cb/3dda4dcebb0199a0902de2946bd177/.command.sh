#!/usr/bin/env Rscript

.libPaths('2')
library(monocle3)

dir.create("temp_fold")
cds <- readRDS("Spleen_cds.RDS")
cell_qc <- read.csv("Spleen_cell_qc.csv")

if(nrow(pData(cds)) > 0) {
    if("false" == 'false') {
        scrublet_out <- read.csv("Spleen_scrublet_out.csv", header=F)
        pData(cds)$scrublet_score <- scrublet_out$V1
        pData(cds)$scrublet_call <- ifelse(scrublet_out$V2 == 1, "Doublet", "Singlet")
        cell_qc$scrublet_score <- scrublet_out$V1
        cell_qc$scrublet_call <- ifelse(scrublet_out$V2 == 1, "Doublet", "Singlet")
    }
}

write.csv(cell_qc, quote=FALSE, file="temp_fold/Spleen_cell_qc.csv")

dup_stats <- read.table(paste0("Spleen", ".duplication_rate_stats.txt"))

df <- data.frame(sample="Spleen", n.reads = dup_stats$V2, n.umi = dup_stats$V3,
                 duplication_rate = dup_stats$V4,
                 doublet_count = sum(cell_qc$scrublet_call == "Doublet", na.rm=TRUE),
                 doublet_perc = paste0(round(sum(cell_qc$scrublet_call == "Doublet",
                                                 na.rm=TRUE)/nrow(cell_qc) * 100, 1), "%"),
                 doublet_NAs=sum(is.na(cell_qc$scrublet_call)))

write.csv(df, file=paste0("Spleen", "_sample_stats.csv"), quote=FALSE, row.names=FALSE)
saveRDS(cds, file="temp_fold/Spleen_cds.RDS")

if ("Spleen" == "Barnyard") {
    fData(cds)$mouse <- grepl("ENSMUSG", fData(cds)$id)
    fData(cds)$human <- grepl("ENSG", fData(cds)$id)

    pData(cds)$mouse_reads <- Matrix::colSums(exprs(cds)[fData(cds)$mouse,])
    pData(cds)$human_reads <- Matrix::colSums(exprs(cds)[fData(cds)$human,])
    pData(cds)$total_reads <- pData(cds)$mouse_reads + pData(cds)$human_reads
    pData(cds)$human_perc <- pData(cds)$human_reads/pData(cds)$total_reads
    pData(cds)$mouse_perc <- pData(cds)$mouse_reads/pData(cds)$total_reads
    pData(cds)$collision <- ifelse(pData(cds)$human_perc >= .9 | pData(cds)$mouse_perc >= .9, FALSE, TRUE)

    collision_rate <- round(sum(pData(cds)$collision/nrow(pData(cds))) * 200, 1)
    fileConn<-file("Barn_collision.txt")
    writeLines(paste0("Spleen", "	", collision_rate, "%"), fileConn)
    close(fileConn)
} else {
    fileConn<-file("Spleen_no_collision.txt")
    writeLines(paste0("Spleen", "	", "NA"), fileConn)
    close(fileConn)
}
