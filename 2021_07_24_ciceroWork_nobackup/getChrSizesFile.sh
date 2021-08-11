#! /bin/bash


# Get the chrom file from the BBI directory used in the ATAC pipeline
cp /net/bbi/vol1/data/genomes_stage/human/human_atac/chromosome_sizes_fasta_finished.txt ./chromosome_sizes_fasta_finished.txt

# Add "chr" to front of the chrom names files
cat ./chromosome_sizes_fasta_finished.txt | awk '$1="chr"$1' | awk '{print $1 "\t" $2}' > hg38chromSizesFormattedForCicero.txt



