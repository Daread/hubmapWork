
library(Matrix)


motifMatFile = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/W145.heart.apex.s1/motif_matrices/W145.heart.apex.s1-peak_motif_matrix.mtx.gz"
motifMat = readMM(gzfile(motifMatFile))

matrixCol = read.table("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/W145.heart.apex.s1/motif_matrices/W145.heart.apex.s1-peak_motif_matrix.columns.txt",
							header=FALSE)
matrixRow = read.table("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/W145.heart.apex.s1/motif_matrices/W145.heart.apex.s1-peak_motif_matrix.rows.txt",
						header=FALSE)

colnames(motifMat) = paste0("chr", matrixCol$V1 )
rownames(motifMat) = matrixRow$V1



# ageLinkedFile = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_10_30_lyonCiceroComp_nobackup/fileInputs/ageDEsites.rds"

ageLinkedFile = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_10_30_lyonCiceroComp_nobackup/fileInputs/sigAgeDEsites.rds"

ageLinkedSites = readRDS(ageLinkedFile)

motifsToCheck = c("NFKB1", "NFKB2")

for (eachMotif in motifsToCheck){

	thisMotifVec = as.vector(motifMat[eachMotif,])
	names(thisMotifVec) = colnames(motifMat)

	# Separate for age-DE-linked or not
	ageDElinkedSites = thisMotifVec[names(thisMotifVec) %in% ageLinkedSites]
	nonDElink = thisMotifVec[!(names(thisMotifVec) %in% ageLinkedSites)]

	ageDElinkedSites = as.numeric(ageDElinkedSites > 0)
	nonDElink = as.numeric(nonDElink > 0)

	# Proportion of DE-linked vs. non with this motif
	propInDELink = sum(ageDElinkedSites) * 1.0 / length(ageDElinkedSites)
	propNonDE    = sum(nonDElink) * 1.0 / length(nonDElink)

	# Output
	print(paste0("Working on ", eachMotif))
	print(paste0("Linked proportion: ", propInDELink))
	print(paste0("Unlinked proportion: ", propNonDE))

	# Test statistic
	myRes = prop.test(c(sum(ageDElinkedSites), sum(nonDElink)), c(length(ageDElinkedSites), length(nonDElink)))
	print(myRes)

}







