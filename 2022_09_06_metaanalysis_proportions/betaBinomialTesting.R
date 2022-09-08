---
title: "Beta binomial cell type abundance testing"
output: html_notebook
date: "9/8/2021"
author: "Lauren Saunders"
---
```{r setup, echo=FALSE, warning=FALSE, message=FALSE}
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_knit$set(progress = TRUE, verbose = FALSE)
suppressPackageStartupMessages({ 
  library(tidyverse)
  library(ggplot2)
  library(data.table)
  library(devtools)
  library(stringr)
  library(reshape2)
  library(VGAM)
  library(glasso)
  library(monocle3)
  library(googledrive)
  library(googlesheets4)
  library(gprofiler2)
  library(gt)
  library(reshape2)
  library(uwot)
  library(gghighlight)
  library(pheatmap)
  library(viridis)
})
```
### Defining the cell types that depend on changes in condition/genotype/treatment/etc.
#### Cellular dependencies
# We can enumerate the cell types that change in frequency in mutant embryos relative to wild-type by modeling cell type abundances with a regression model.
# ```{r bb help functions}
# a help function to tidy the vgam model output - used in compare_abundance function
tidy.vglm = function(x, conf.int=FALSE, conf.level=0.95) {
  co <- as.data.frame(coef(summary(x)))
  names(co) <- c("estimate","std.error","statistic","p.value")
  if (conf.int) {
    qq <- qnorm((1+conf.level)/2)
    co <- transform(co,
                    conf.low=estimate-qq*std.error,
                    conf.high=estimate+qq*std.error)
  }
  co <- data.frame(term=rownames(co),co)
  rownames(co) <- NULL
  return(co)
}
# function that uses beta binomial regression to compare the abundances of a celltype from two conditions
compare_abundance = function(cell_df, wt_df, model_formula = "count_df ~ genotype", ...){
  comb_df = rbind(as.data.frame(cell_df), as.data.frame(wt_df)) %>% 
    mutate(genotype = fct_relevel(genotype, c(unique(wt_df$genotype))))
  
  print(unique(comb_df$genotype))
  
  cell.types = unique(comb_df$cell_type)
  
  test_res = sapply(cell.types, 
                    FUN = function(x) {
                      type_df = comb_df %>% filter(cell_type == x)
                      count_df = cbind(type_df$cells, type_df$total_cells - type_df$cells)
                      print(head(type_df))
                      fit =  vglm(as.formula(model_formula), betabinomial, data = type_df, trace = TRUE, ...)
                      fit_df = tidy.vglm(fit)}, USE.NAMES = T, simplify = F)
  
  test_res = do.call(rbind, test_res)
  test_res = test_res %>% tibble::rownames_to_column(var = "cell_group")
  test_res %>% arrange(desc(estimate))
}
# ```
# First we will generate cds objects for cell abundances
# ```{r load data}
# load coldata with cell type annotations
coldata_df = fread("/Volumes/GoogleDrive/Shared drives/Trapnell Lab/Projects/SDG/GAPFISH/data/final_final_ref/R_objects_other/gap16_all-cells_anno_clean_2.7M_coldata.csv",
                       sep = ",", stringsAsFactors = F, data.table = F, na.strings = "") %>% 
  mutate(timepoint = as.character(timepoint))
# timepoint vector
timepoints = sort(unique(coldata_df$timepoint))
timepoints
# ```
# ```{r compute size factors by timepoint}
cds_list = list()
for (time in timepoints){
  message(paste0("making cell count CDS for ", time))
  
  # generate counts for all 18h fish
  coldata_sub = coldata_df %>% 
    filter(timepoint == time)
  
  counts_df = coldata_sub %>% 
    group_by(embryo, cell_type_broad) %>%
    filter(cell_type_broad != "") %>%
    dplyr::summarize(cells=n()) %>%
    ungroup() %>% 
    pivot_wider(names_from = cell_type_broad, values_from = cells, values_fill = c(0))
  
  count_mat = as.matrix(counts_df[,-1])
  rownames(count_mat) = counts_df$embryo
  
  # grab meta data
  meta_df = coldata_sub %>%
    select(embryo, timepoint, gene_target) %>% 
    distinct() %>% 
    ungroup()
  
  # counts
  rownames(meta_df) = meta_df$embryo
  count_mat = count_mat[as.character(meta_df$embryo),] # same order
  
  cell_cds = new_cell_data_set(t(count_mat), 
                               cell_metadata=meta_df) %>% 
    preprocess_cds(num_dim = 10, 
                   norm_method="size_only", 
                   method = "PCA")
  
  cds_list[[time]] <- cell_cds
}
# get size factors per emb for normalization in next step
groups = names(cds_list) %>% sort()
sf_list = list()
for (i in groups) {
  sf_df = size_factors(cds_list[[i]]) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column()
  colnames(sf_df) = c("embryo" ,"size_factor")
  
  sf_list[[i]] = sf_df
}
size_factor_df = do.call(rbind, sf_list)
rownames(size_factor_df) <- NULL
fwrite(size_factor_df, "all-gap16_update_per-timepoint_size-factor_df.csv", sep = ",")
```
```{r summarize cell counts across embryos}
# if already computed, load size factor df
size_factor_df = fread("/Volumes/GoogleDrive/Shared drives/Trapnell Lab/Projects/SDG/GAPFISH/data/final_final_ref/pre_comp/all-gap16_update-emb_per-timepoint_size-factor_df.csv",
                       sep = ",", stringsAsFactors = F, data.table = F)
# generate coldata list by timepoint
coldat_list = lapply(timepoints,
                     function(x) {
                       sub_coldat = coldata_df %>% 
                         filter(timepoint == x)
                       sub_coldat
                     })
names(coldat_list) <- timepoints
# generate normalized summary tables
annos = "cell_type_broad"
summary_list = list()
for (time in timepoints) {
  message(paste0("making cell summary for ", time, "hpf"))
  
  sub_df = coldat_list[[time]]
  all_summary = sub_df %>%
    group_by(genotype = gene_target, embryo, cell_type = get(annos)) %>% 
    dplyr::summarize(cells = n(), .groups = "keep") %>%
    left_join(size_factor_df, by = "embryo") %>% 
    mutate(cells = round(cells / size_factor)) %>% 
    select(genotype, embryo, cell_type, cells) %>% 
    ungroup() %>%
    distinct() %>%
    group_by(embryo) %>% 
    mutate(total_cells = sum(cells), 
           cell_type = as.character(cell_type)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = cell_type, values_from = cells, values_fill = list(cells = 0)) %>%
    pivot_longer(cols = -c(genotype, embryo, total_cells), values_to = "cells", names_to = "cell_type")
  summary_list[[time]] <- all_summary
}
# collapse and save summary list
summary_merge = rbindlist(summary_list, idcol = "timepoint") %>% 
  as.data.frame()
fwrite(summary_merge, "saved_celltype_summary_list.csv", sep = ",", na = "NA")
```
```{r compute embryo means and fold changes}
# make cell type means/embryo lists
celltype_means_list = list()
for (time in timepoints){
  mean_df = summary_list[[time]] %>% 
    group_by(cell_type) %>%
    summarize(emb_mean = round(mean(cells))) %>%
    ungroup()
  
  celltype_means_list[[time]] = mean_df
}
# compute log2fc relative to ctrl for all celltypes for each timepoint
fc_summary_list <- list()
for (time in timepoints){
  
  # filter geno groups
  gene_targets <- unique(coldat_list[[time]]$gene_target)
  ctrl_groups <- grep("ctrl", gene_targets, value = TRUE)
  test_groups <- grep("ctrl", gene_targets, invert = TRUE, value = TRUE)
  
  # get cell summary
  all_summary = summary_list[[time]]
  
  # wildtype summary table
  wt_df = all_summary %>%
    filter(genotype == "ctrl-inj") %>%
    group_by(genotype, cell_type) %>% 
    summarize(ctrl_mean = round(mean(cells))) %>% 
    ungroup() %>% 
    distinct(cell_type, ctrl_mean)
  
  # perturbation summary table
  fc_df = all_summary %>%
    filter(genotype != "ctrl-inj") %>% 
    group_by(genotype, cell_type) %>% 
    summarize(geno_mean = round(mean(cells))) %>% 
    left_join(wt_df, by = "cell_type") %>% 
    mutate(abund_log2fc = log2((geno_mean + 1)/(ctrl_mean+1)))
  
  fc_summary_list[[time]] <- fc_df
}
```
```{r run BB testing across cell types and timepoints}
mean_thresh = 4 # threshold for mean number of cells per embryo for whether to test a cell type
res_list <- list()
for (time in timepoints){
  
  # filter geno groups
  gene_targets <- unique(coldat_list[[time]]$gene_target)
  
  # get cell summary
  all_summary = summary_list[[time]]
  
  # filter cell groups
  cell.groups = all_summary %>% 
    group_by(cell_type) %>% 
    summarize(emb_mean = mean(cells)) %>%
    ungroup() %>% 
    filter(emb_mean > mean_thresh) %>% 
    pull(cell_type)
  
  cell.groups <- as.character(cell.groups)
  
  message(paste0("comparing ", length(cell.groups), " cell types at ", time))
  
  # filter for control groups only
  wt_filt = all_summary %>%
    dplyr::filter(genotype == "ctrl-inj") %>%
    dplyr::filter(cell_type %in% cell.groups)
  
  # filter for perturbation groups only
  mut_filt = all_summary %>%
    filter(genotype != "ctrl-inj") %>%
    dplyr::filter(cell_type %in% cell.groups)
  
  celltype_diff = compare_abundance(mut_filt, wt_filt, model_formula = "count_df ~ genotype") %>% 
    filter(!(grepl("intercept", term, ignore.case = T))) %>% 
    dplyr::mutate(qval = p.adjust(p.value, method = "BH")) %>% 
    mutate(geno.v.wt = stringr::str_sub(term, 9)) %>% 
    separate(cell_group, into = c("cell_group", NULL), sep = "\\.") %>%
    select(geno.v.wt, everything(), -term) %>%
    arrange(geno.v.wt) %>% 
    left_join(celltype_means_list[[time]] %>% 
                dplyr::rename(cell_group = cell_type), by = c("cell_group")) %>%
    left_join(fc_summary_list[[time]] %>% 
                select(geno.v.wt = genotype, cell_group = cell_type, abund_log2fc),
              by = c("geno.v.wt", "cell_group"))
  
  celltype_diff$timepoint = time
  
  res_list[[time]] <- celltype_diff
}
# save res list
saveRDS(res_list, "bb_results_res_list.RDS")
# collapse list into dataframe
res_df = do.call(rbind, res_list_peri)
# save results
fwrite(res_df, "bb_res_df_date.csv", sep = ",", na = "NA")
```
```{r plot heatmaps across timepoints}
# set qval threshold
qval_thresh = 0.05
anno_type = "broadtype"
dir_name = "per_timepoint_heatmaps"
dir.create(dir_name) # make a new directory
for (time in timepoints){
  emb_df = celltype_means_list[[time]]
  emb_df = as.data.frame(emb_df)
  rownames(emb_df) = emb_df$cell_type
  emb_df = emb_df %>%
    select(emb_mean)
  hm_df = res_list[[time]] %>% 
    ungroup() %>% 
    select(gene_target = geno.v.wt, cell_group, abund_log2fc, qval) %>%
    dplyr::mutate(sig_fc = case_when(qval < qval_thresh ~ abund_log2fc,
                                      TRUE ~ 0)) %>%
    select(gene_target, cell_group, sig_fc) %>%
    pivot_wider(names_from = cell_group,
                values_from = sig_fc)
  
  hm_mat = as.matrix(hm_df[,-1])
  rownames(hm_mat) = hm_df$gene_target
  hm_mat[is.na(hm_mat)] <- 0
  hm_mat = hm_mat[,Matrix::colSums(abs(hm_mat)) > 0]
  
  # even out boundaries
  cutoff = floor(2*min(abs(min(hm_mat)), max(hm_mat))) / 2 # round to nearest half step down
  hm_mat[hm_mat > cutoff] <- cutoff
  hm_mat[hm_mat < -cutoff] <- -cutoff
  
  print(cutoff)
  message(paste0("plotting heatmap for abundance changes in ", time, "hpf"))
  
  pheatmap::pheatmap(hm_mat, cluster_cols = T, cluster_rows = T,
                     cellheight = 10, cellwidth = 10, fontsize=8, border_color = NA,
                     annotation_col = log10(emb_df), #annotation_names_col = F,
                     color = colorRampPalette(colors = c("darkblue", "white", "red4"))(100),
                     angle_col = 45,
                     filename = paste0(dir_name,"/", time, "_all-geno_", anno_type, "_heatmap_update.png"))
}
```
```{r plot boxplots of abundance changes for specific types and perturbations}
# make new directory
dir_name = "celltype_boxplots"
dir.create(dir_name)
# plot one cell type over time
types = c("notochord")
genos = c("ctrl-inj", "tbxta", "noto")
times = c("18", "24", "36")
summary_merge %>% 
  filter(genotype %in% genos) %>%
  filter(cell_group %in% types) %>% 
  filter(timepoint %in% times) %>%
  ggplot(aes(x = factor(genotype, levels = genos), 
             y = cells, 
             fill = genotype)) +
  geom_boxplot(color = "black", size = 0.3, outlier.size = 0.5, outlier.shape = NA) +
  geom_jitter(size = 0.5, color = "black") +
  theme(axis.title.x = element_blank(), legend.position = "none", 
        axis.ticks = element_line(color = "black", size = .3),
        legend.title = element_blank(), 
        axis.text.x = element_text(angle= 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black")) +
  facet_wrap(~cell_type+timepoint, nrow = 2, scales = "free_y") +
  monocle3:::monocle_theme_opts()
ggsave(paste0(dir_name, paste0(types, "_ctrl-tbxta-noto_abund_boxplot.png")), dpi = 500, width = 6, height = 3)
```
