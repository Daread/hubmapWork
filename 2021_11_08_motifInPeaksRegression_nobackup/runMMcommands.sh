


# 11-8-21: Run models using age on a linear scale. Using a log scale immediately raise eyebrows when I spoke w/ Jennifer, 
#          and probably will be a thorn in my side to justify to reviewers. If they ask to compare likelihoods/errors I can, but
#          probably simplest to just analyze based on linear age

qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Endocardium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Fibroblast Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Macrophage Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Perivascular Anatomical_Site,Age,Sex





# 9-16-21: Running code again but with a tweak to 
# qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Age,Sex
# qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Endocardium Anatomical_Site,Age,Sex
# qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Fibroblast Anatomical_Site,Age,Sex
# qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Macrophage Anatomical_Site,Age,Sex
# qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh T_Cell Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Adipocytes Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh B_Cell Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Mast_Cell Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Neuronal Anatomical_Site,Age,Sex
# qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh VSM_and_Pericyte Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Lymphatic_Endothelium Anatomical_Site,Age,Sex



# # 11-22-21 Plotting
# Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult

# Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript plotOutputs.R -c Perivascular -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript plotOutputs.R -c Macrophage -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult

# Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult



# # 9-20-21 Plotting
# Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript plotOutputs.R -c Macrophage -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript plotOutputs.R -c T_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript plotOutputs.R -c VSM_and_Pericyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult



# # 9-21-21 Add GSEA work
# Rscript gseaAnalysis.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript gseaAnalysis.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript gseaAnalysis.R -c Macrophage -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript gseaAnalysis.R -c T_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript gseaAnalysis.R -c VSM_and_Pericyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript gseaAnalysis.R -c Fibroblast -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult


# # 9-21-21 Add GSEA work
# Rscript gseaAnalysis.R -c Endocardium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript gseaAnalysis.R -c B_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript gseaAnalysis.R -c Mast_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript gseaAnalysis.R -c Neuronal -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript gseaAnalysis.R -c Lymphatic_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult
# Rscript gseaAnalysis.R -c Adipocytes -m _fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult










# 12-3-21: Run regression by cell type to find cell-type spec motifs
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Endocardium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Fibroblast Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Macrophage Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh VSM_and_Pericyte Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh T_Cell Anatomical_Site,Age,Sex



# 12-3-21: Run regression within cell types to see what motifs change by age/sex in those cells
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Endocardium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Fibroblast Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh Macrophage Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh VSM_and_Pericyte Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixed.sh T_Cell Anatomical_Site,Age,Sex


# 12-8-21 Get the outputs

Rscript plotOutputs.R -c T_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult
Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult
Rscript plotOutputs.R -c Macrophage -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult
Rscript plotOutputs.R -c VSM_and_Pericyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult
Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult
Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult


# Cell type specific work
Rscript plotOutputs.R -c Cell_Type_Specificity_forT_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/
Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Macrophage -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c VSM_and_Pericyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for





# 1-21-21: 
# I've re-run the code that generates peak x motif matrices to use different p values as the FIMO motif cutoff. 
# Modified code to re-run and test the results:

qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixedAndP.sh Cardiomyocyte Anatomical_Site,Age,Sex 1e-6
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixedAndP.sh Vascular_Endothelium Anatomical_Site,Age,Sex 1e-6
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixedAndP.sh Endocardium Anatomical_Site,Age,Sex 1e-6
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixedAndP.sh Fibroblast Anatomical_Site,Age,Sex 1e-6
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixedAndP.sh Macrophage Anatomical_Site,Age,Sex 1e-6
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixedAndP.sh VSM_and_Pericyte Anatomical_Site,Age,Sex 1e-6
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleMMrunSetFixedAndP.sh T_Cell Anatomical_Site,Age,Sex 1e-6



# 1-23-21: See the output results
Rscript plotOutputs.R -c T_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_p1e-6_MMresult
Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_p1e-6_MMresult
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_p1e-6_MMresult
Rscript plotOutputs.R -c Macrophage -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_p1e-6_MMresult
Rscript plotOutputs.R -c VSM_and_Pericyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_p1e-6_MMresult
Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_p1e-6_MMresult
Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_p1e-6_MMresult










# 2-7-22: Run regression by cell type to find cell-type spec motifs, adding UMI to the 
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Age,Sex,log10umi
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Age,Sex,log10umi
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Endocardium Anatomical_Site,Age,Sex,log10umi
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Fibroblast Anatomical_Site,Age,Sex,log10umi
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Macrophage Anatomical_Site,Age,Sex,log10umi
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh VSM_and_Pericyte Anatomical_Site,Age,Sex,log10umi
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh T_Cell Anatomical_Site,Age,Sex,log10umi




# Cell type specific work re-run plot code to get csv with all test results
# Rscript plotOutputs.R -c Cell_Type_Specificity_forT_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/
Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Macrophage -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c VSM_and_Pericyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for


# 2-7-22: Run regression by cell type to find cell-type spec motifs, adding UMI and FRIT to the covariates
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Age,Sex,log10umi,FRIT
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Age,Sex,log10umi,FRIT
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Endocardium Anatomical_Site,Age,Sex,log10umi,FRIT
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Fibroblast Anatomical_Site,Age,Sex,log10umi,FRIT
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Macrophage Anatomical_Site,Age,Sex,log10umi,FRIT
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh VSM_and_Pericyte Anatomical_Site,Age,Sex,log10umi,FRIT
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh T_Cell Anatomical_Site,Age,Sex,log10umi,FRIT


# 2-7-22: Run regression by cell type. ONLY using UMI as a fixed effect covariate (in addition to being the cell type of interest/not). Trying to exactly mirror the fetal atlas setup
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Cardiomyocyte log10umi
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Vascular_Endothelium log10umi
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Endocardium log10umi
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Fibroblast log10umi
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh Macrophage log10umi
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh VSM_and_Pericyte log10umi
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixed.sh T_Cell log10umi




# ,log10umi

# Cell type specific work, 2-7-22 adding in log10umi
Rscript plotOutputs.R -c T_Cell -m _fix_Anatomical_Site,Age,Sex,log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age,Sex,log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age,Sex,log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Macrophage -m _fix_Anatomical_Site,Age,Sex,log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c VSM_and_Pericyte -m _fix_Anatomical_Site,Age,Sex,log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Age,Sex,log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Age,Sex,log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for



# Cell type specific work, 2-7-22 adding in log10umi and FRIT
Rscript plotOutputs.R -c T_Cell -m _fix_Anatomical_Site,Age,Sex,log10umi,FRIT_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age,Sex,log10umi,FRIT_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age,Sex,log10umi,FRIT_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Macrophage -m _fix_Anatomical_Site,Age,Sex,log10umi,FRIT_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c VSM_and_Pericyte -m _fix_Anatomical_Site,Age,Sex,log10umi,FRIT_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Age,Sex,log10umi,FRIT_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Age,Sex,log10umi,FRIT_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for




# Cell type specific work, 2-7-22 adding in log10umi and FRIT
Rscript plotOutputs.R -c T_Cell -m _fix_log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Macrophage -m _fix_log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c VSM_and_Pericyte -m _fix_log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Fibroblast -m _fix_log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Endocardium -m _fix_log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult -p cellTypeSpec/Cell_Type_Specificity_for






# 2-7-22: Run regression by cell type to find cell-type spec motifs, adding UMI to the 
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixedPeakTypes.sh Cardiomyocyte Anatomical_Site,Age,Sex,log10umi
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixedPeakTypes.sh Vascular_Endothelium Anatomical_Site,Age,Sex,log10umi
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixedPeakTypes.sh Endocardium Anatomical_Site,Age,Sex,log10umi
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixedPeakTypes.sh Fibroblast Anatomical_Site,Age,Sex,log10umi
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixedPeakTypes.sh Macrophage Anatomical_Site,Age,Sex,log10umi
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/singleCellMMrunSetFixedPeakTypes.sh Perivascular Anatomical_Site,Age,Sex,log10umi







