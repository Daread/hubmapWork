



# 2-7-22: Run regression by cell type to find cell-type spec motifs, adding UMI to the 
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh Endocardium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh Fibroblast Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh Macrophage Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh VSM_and_Pericyte Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh T_Cell Anatomical_Site,Age,Sex


qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh Endocardium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh Fibroblast Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh Macrophage Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh VSM_and_Pericyte Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh T_Cell Anatomical_Site,Age,Sex


# 3-15-22: Run regression with smaller cell types
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh Adipocytes Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh B_Cell Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh Lymphatic_Endothelium Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh Mast_Cell Anatomical_Site,Age,Sex
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleCellMMrunSetFixed.sh Neuronal Anatomical_Site,Age,Sex




qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh Adipocytes Anatomical_Site,Age,Sex 
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh B_Cell Anatomical_Site,Age,Sex 
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh Lymphatic_Endothelium Anatomical_Site,Age,Sex 
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh Mast_Cell Anatomical_Site,Age,Sex 
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/singleMMrunSetFixed.sh Neuronal Anatomical_Site,Age,Sex 



# Cell_Type_Specificity_forEndocardium_fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult
# Fibroblast_fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult


# Cell type specific work, 2-7-22 adding in log10umi
Rscript plotOutputs.R -c Macrophage -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c T_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c VSM_and_Pericyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for

# 3-16-22: Plot the outputs of the cell-specific testing for rarer types
Rscript plotOutputs.R -c Adipocytes -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c B_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Lymphatic_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Mast_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Neuronal -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for




# Cell type specific work, 2-7-22 adding in log10umi
Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult  
Rscript plotOutputs.R -c Macrophage -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult  
Rscript plotOutputs.R -c T_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult  
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult  
Rscript plotOutputs.R -c VSM_and_Pericyte -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult  
Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult  
Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult  


# Adding in regression 





# 3-16-22: Plot the outputs of the cell-specific testing for rarer types
Rscript plotOutputs.R -c Adipocytes -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c B_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Lymphatic_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Mast_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for
Rscript plotOutputs.R -c Neuronal -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult -p cellTypeSpec/Cell_Type_Specificity_for


# 



# 4-22-22: Plot the outputs of the testing for rarer types
Rscript plotOutputs.R -c Adipocytes -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult 
Rscript plotOutputs.R -c B_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult 
Rscript plotOutputs.R -c Lymphatic_Endothelium -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult
Rscript plotOutputs.R -c Mast_Cell -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult
Rscript plotOutputs.R -c Neuronal -m _fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult 






