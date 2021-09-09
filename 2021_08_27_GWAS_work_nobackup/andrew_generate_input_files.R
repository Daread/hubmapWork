library(stringr)
library(readr)
library(dplyr)
library(glue)
options(stringsAsFactors=FALSE)
args = list()
#args$specificity_scores = "/net/trapnell/vol1/ajh24/proj/2018three_level_sciatac/results/ahill/2019_05_06_snapatac_exploratory_analysis/differential_tests/all_tissues/peak_specificity_scores.annotation.txt"
#args$specificity_scores = "/net/trapnell/vol1/ajh24/proj/2018three_level_sciatac/results/ahill/2019_05_06_snapatac_exploratory_analysis/differential_tests/all_tissues.no_tissue_label/peak_specificity_scores.annotation.txt"
#args$specificity_scores = "/net/shendure/vol10/projects/silvia/BDRL_atac/analysis/clustering/nobackup/nnls/annotation_2zwindows/peak_specificity_scores.54celltypes.txt"
#args$specificity_scores = "/net/trapnell/vol1/ajh24/proj/2018three_level_sciatac/results/ahill/2019_05_06_snapatac_exploratory_analysis/differential_tests/brain/peak_specificity_scores.cluster.txt"
args$specificity_scores = "/net/trapnell/vol1/ajh24/proj/2018three_level_sciatac/results/ahill/2019_05_06_snapatac_exploratory_analysis/differential_tests/all_tissues/peak_specificity_scores.cluster.txt"
args$top_peaks = 10000

#args$output_dir = 'top_peak_bed_files.final_annotations'
#args$sample_sheet_out = 'samplesheet.final_annotations.txt'
#args$output_dir = 'top_peak_bed_files.brain_clusters'
args$output_dir = "top_peak_bed_files.fetal_specific"

args$sample_sheet_out = 'samplesheet.fetal_specific.txt'
args$peak_metadata = '/net/shendure/vol10/projects/silvia/BDRL_atac/analysis/clustering/nobackup/nnls/annotation_2zwindows/peak_metadata.txt'
args$peak_whitelist = 'fetal_specific_peaks.txt'
peak_whitelist = readr::read_delim(args$peak_whitelist, delim=',', col_names=c('peak'))

#specificity_scores = readr::read_delim('/net/trapnell/vol1/ajh24/proj/2018three_level_sciatac/results/ahill/2019_05_06_snapatac_exploratory_analysis/differential_tests/all_tissues/peak_specificity_scores.annotation.txt', '\t')
specificity_scores = readr::read_delim(args$specificity_scores, delim='\t')

message('Subsetting to overlapping encode fetal tissues...')
specificity_scores = subset(specificity_scores, grepl('adrenal|spleen|stomach|intestine|heart', group))

peak_metadata = read.table(args$peak_metadata, sep='\t')
peak_metadata = subset(peak_metadata,!zscore_filtered)
specificity_scores = subset(specificity_scores, feature %in% peak_metadata$peak) 
specificity_scores = subset(specificity_scores, feature %in% peak_whitelist$peak)

# Now get the top peaks per annotation
specificity_scores = specificity_scores %>%
    group_by(group) %>%
    filter(specificity_score > 0) %>% # TODO think this is right choice
    arrange(desc(specificity_score)) %>%
    slice(1:args$top_peaks) %>%
    ungroup()

specificity_scores$group = str_replace_all(specificity_scores$group, ' ', '_')
specificity_scores$group = str_replace_all(specificity_scores$group, '-', '_')
specificity_scores$group = str_replace_all(specificity_scores$group, '\\?', '_unsure')
specificity_scores$group = str_replace_all(specificity_scores$group, '/', '_')

# Write each out to a file
sample_sheet_rows = lapply (unique(specificity_scores$group), function(annotation) {
    output_file = glue('{args$output_dir}/{annotation}.bed')
    data.frame(sample_id=as.character(annotation), sites=as.character(output_file))
})

sample_sheet_df = bind_rows(sample_sheet_rows)
readr::write_delim(sample_sheet_df, path=args$sample_sheet_out, delim='\t')

dir.create(args$output_dir, showWarnings=FALSE)
for (entries in sample_sheet_rows) {
    annotation = entries$sample_id[1]
    output_file = entries$sites[1]

    peaks = filter(specificity_scores, group == annotation)
    location = str_split_fixed(peaks$feature, '_', 3)
    output_dataframe = data.frame(chrom=location[, 1], start=as.numeric(location[, 2]), stop=as.numeric(location[, 3]))
    output_dataframe = arrange(output_dataframe, chrom, start)
    readr::write_delim(output_dataframe, delim='\t', path=output_file, col_names=FALSE)
}
