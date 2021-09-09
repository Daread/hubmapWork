library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(viridis)
library(tidyr)
library(ggrepel)
options(stringsAsFactors=FALSE)

##########################################
# Heatmap
##########################################
# Threshold for the final heatmap (which points get shown)
FINAL_QVAL_THRESHOLD = 0.2
FINAL_ENRICHMENT_THRESHOLD = 0

# These are just two different sets of traits

input_file = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/fileOutputs/bedFilesForGWAS_1000_200_0.05_5_600/ldsc_results/results_gathered.txt"
# input_file = "/net/trapnell/vol1/ajh24/proj/2018three_level_sciatac/results/ahill/2019_10_02_initial_gwas_integration/ldsc_results.ukbiobank_alkes_price_409k.final_annotations/results_gathered.txt"
#input_file = "/net/trapnell/vol1/ajh24/proj/2018three_level_sciatac/results/ahill/2019_10_02_initial_gwas_integration/ldsc_results.final_annotations/results_gathered.txt"

# Read in results and modify the labels from file to make shorter
ukbb_results = readr::read_delim(input_file, delim='\t') %>%
                mutate(result_file = str_replace(result_file, 'ldsc_results.[^/]+/score_sumstats/', '')) %>%
                mutate(result_file = str_replace(result_file, '.results', '')) %>%
                mutate(result_file = str_replace(result_file, 'blood_', '')) %>%
                mutate(result_file = str_replace(result_file, 'disease_', '')) %>%
                mutate(result_file = str_replace(result_file, 'repro_', '')) %>%
                mutate(result_file = str_replace(result_file, 'bp_', '')) %>%
                mutate(result_file = str_replace(result_file, 'lung_', '')) %>%
                mutate(result_file = str_replace(result_file, 'other_', '')) %>%
                mutate(result_file = str_replace(result_file, 'mental_', '')) %>%
                mutate(result_file = str_replace(result_file, 'cov_', '')) %>%
                mutate(result_file = str_replace(result_file, 'body_', '')) %>%
                mutate(result_file = str_replace(result_file, 'impedence_', '')) %>%
                mutate(result_file = str_replace_all(result_file, '_', ' ')) %>%
                mutate(cell_type = str_split_fixed(result_file, '-', n=2)[,1]) %>%
                mutate(cell_type = str_replace_all(cell_type, '_', ' ')) %>%
                mutate(trait = str_split_fixed(result_file, '-', n=2)[,2]) %>%
                dplyr::mutate(qval = p.adjust(Enrichment_p, method='BH'))

# Only consider traits with at least one significant positive enrichment and reasonable heritability
traits_to_consider = ukbb_results %>%
                    dplyr::filter(qval < FINAL_QVAL_THRESHOLD & Enrichment > FINAL_ENRICHMENT_THRESHOLD & h2 > 0.01)

ukbb_results.filtered = dplyr::filter(ukbb_results, trait %in% traits_to_consider$trait)
significant_plot_df = subset(ukbb_results.filtered, qval < FINAL_QVAL_THRESHOLD & Enrichment > FINAL_ENRICHMENT_THRESHOLD)

significant_plot_df = significant_plot_df[!duplicated(significant_plot_df %>% select(trait, cell_type)),]
gwas_cluster_matrix = as.data.frame(spread(significant_plot_df %>% select(trait, cell_type, Enrichment), "cell_type", "Enrichment"))
row.names(gwas_cluster_matrix) = gwas_cluster_matrix[, 1]
gwas_cluster_matrix = gwas_cluster_matrix[, -1]
gwas_cluster_matrix[is.na(gwas_cluster_matrix)] = 0
hc_cols_results = hclust(dist(gwas_cluster_matrix), method='ward.D2')
hc_rows_results = hclust(dist(t(gwas_cluster_matrix)), method='ward.D2')

significant_plot_df$trait = factor(significant_plot_df$trait, levels = rownames(gwas_cluster_matrix)[hc_cols_results$order])
significant_plot_df$cell_type = factor(significant_plot_df$cell_type, levels = colnames(gwas_cluster_matrix)[hc_rows_results$order])

# One option is the heatmap
ggplot(significant_plot_df, aes(cell_type, trait, fill=Enrichment)) +
    geom_tile(color='black', fill='white') +
    geom_tile(color='black', aes(fill=log10(Enrichment))) +
    theme_classic() +
    scale_fill_viridis(option='magma', name='log10(Enrichment)') +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    xlab('trait') +
    ylab('annotation') +
    ggsave('trait_heatmap.png', height=5, width=6.5)

# 10 by 8 works for first trait set (UKBB)
# 5 by 6.5 works for the other trait set

##########################################
# Example of individual volcano plot
##########################################
ukbb_results.single_trait = ukbb_results %>%
    filter(trait == 'MENOPAUSE AGE')

LABEL_QVAL_THRESHOLD = 0.50
LABEL_ENRICHMENT_THRESHOLD = 0

ggplot(ukbb_results.single_trait, aes(Enrichment, -log10(qval))) +
    geom_point() +
    geom_hline(yintercept=-log10(FINAL_QVAL_THRESHOLD), color='red', linetype='dashed') +
    geom_label_repel(data=subset(ukbb_results.single_trait, qval< LABEL_QVAL_THRESHOLD & Enrichment > LABEL_ENRICHMENT_THRESHOLD), aes(label=cell_type)) +
    theme_classic() +
    ggsave('trait_heatmap.png')
