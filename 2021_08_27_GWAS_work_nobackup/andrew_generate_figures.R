
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(viridis)
library(tidyr)
library(ggrepel)
options(stringsAsFactors=FALSE)

# 9-13-21 David: Adding command line parser so I can specify the hyperparameter set used and use that to define input files/output names
library(optparse)
# Get the passed parameters
option_list = list(
  make_option(c("-u", "--promoterUpstream"), type="numeric", 
        default=0,
              help="Bases upstream of TSS to use for input features", metavar="numeric"),
  make_option(c("-d", "--promoterDownstream"), type="numeric", 
        default=0,
              help="Bases downstream of TSS to use for input features", metavar="numeric"),
  make_option(c("-k", "--coaccessCutoff"), type="numeric", 
        default=.05,
              help="Cutoff for keeping coaccessible peaks", metavar="numeric"),
  make_option(c("-n", "--maxNdistalSites"), type="numeric", 
        default=20,
              help="Max number of sites to link to a gene's promoter", metavar="numeric"),
  make_option(c("-q", "--peakSize"), type="numeric", 
        default=600,
              help="-1 -> Keep as is, or enter a positive integer to make all peaks the same size", metavar="numeric"),

  make_option(c("-m", "--cpmMin"), type="numeric", 
        default=2,
              help="Min log2 CPM to count as expressed in a cell type", metavar="numeric"),

  make_option(c("-r", "--ratMin"), type="numeric", 
        default=.1,
              help="Minimum log2(cpm/averageCPM) value for genes to count in a cell type", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

opt$variableParams = paste0(as.character(opt$promoterUpstream), "_", as.character(opt$promoterDownstream), "_",
                            as.character(opt$coaccessCutoff), "_", as.character(opt$maxNdistalSites), "_",
                            as.character(opt$peakSize), "_",
                             as.character(opt$cpmMin), "_", as.character(opt$ratMin))



##########################################
# Heatmap
##########################################
# Threshold for the final heatmap (which points get shown)
FINAL_QVAL_THRESHOLD = 0.05
FINAL_ENRICHMENT_THRESHOLD = 0

# These are just two different sets of traits
gwasRun = "ldsc_results_full"
input_file = paste0("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/fileOutputs/bedFilesForGWAS_",
                     opt$variableParams, "/", 
                    gwasRun, "/results_gathered.txt")

# Read in results and modify the labels from file to make shorter
outFileLabel = paste0("./fileOutputs/bedFilesForGWAS_", opt$variableParams, "/")

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
                mutate(result_file = str_replace(result_file, outFileLabel, '')) %>% # DFR Add
                mutate(result_file = str_replace(result_file, 'HEIGHTz', 'HEIGHT')) %>% # DFR Add
                mutate(result_file = str_replace(result_file, 'adjMEDz', ' BP')) %>% # DFR Add
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
png(paste0('./plots/trait_heatmap_allCellTypes', opt$variableParams, "q_", as.character(FINAL_QVAL_THRESHOLD), '.png'), height=1000, width=1400, res=200)
myPlot = ggplot(significant_plot_df, aes(cell_type, trait, fill=Enrichment)) +
    geom_tile(color='black', fill='white') +
    geom_tile(color='black', aes(fill=log10(Enrichment))) +
    theme_classic() +
    scale_fill_viridis(option='magma', name='log10(Enrichment)') +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) #+
    # xlab('trait') +
    # ylab('annotation')# 
print(myPlot)
dev.off()

# Explicitly plot q values, as in the mouse atac atlas
png(paste0('./plots/trait_pval_heatmap_allCellTypes', opt$variableParams, "q_", as.character(FINAL_QVAL_THRESHOLD), '.png'), height=1000, width=1400, res=200)
myPlot = ggplot(significant_plot_df, aes(cell_type, trait, fill=Enrichment_p)) +
    geom_tile(color='black', fill='white') +
    geom_tile(color='black', aes(fill=-log10(Enrichment_p))) +
    theme_classic() +
    scale_fill_viridis(option='magma', name='-log10(P_Value)') +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) #+
    # xlab('trait') +
    # ylab('annotation')# 
print(myPlot)
dev.off()



# # One option is the heatmap
# png(paste0('./plots/trait_heritability_heatmap_allCellTypes', opt$variableParams, "q_", as.character(FINAL_QVAL_THRESHOLD), '.png'), height=1000, width=1400, res=200)
# myPlot = ggplot(significant_plot_df, aes(cell_type, trait, fill=Prop._h2)) +
#     geom_tile(color='black', fill='white') +
#     geom_tile(color='black', aes(fill=log10(Prop._h2))) +
#     theme_classic() +
#     scale_fill_viridis(option='magma', name='log10(Proportion_H2)') +
#     coord_flip() +
#     theme(axis.text.x = element_text(angle = 60, hjust = 1)) #+
#     # xlab('trait') +
#     # ylab('annotation')# 
# print(myPlot)
# dev.off()



# 10 by 8 works for first trait set (UKBB)
# 5 by 6.5 works for the other trait set

# ##########################################
# # Example of individual volcano plot
# ##########################################
# ukbb_results.single_trait = ukbb_results %>%
#     filter(trait == 'MENOPAUSE AGE')

# LABEL_QVAL_THRESHOLD = 0.50
# LABEL_ENRICHMENT_THRESHOLD = 0

# ggplot(ukbb_results.single_trait, aes(Enrichment, -log10(qval))) +
#     geom_point() +
#     geom_hline(yintercept=-log10(FINAL_QVAL_THRESHOLD), color='red', linetype='dashed') +
#     geom_label_repel(data=subset(ukbb_results.single_trait, qval< LABEL_QVAL_THRESHOLD & Enrichment > LABEL_ENRICHMENT_THRESHOLD), aes(label=cell_type)) +
#     theme_classic() +
#     ggsave('trait_heatmap.png')
