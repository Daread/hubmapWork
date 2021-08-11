

qsub -P trapnelllab -l mfree=10G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrun.sh Endocardium
qsub -P trapnelllab -l mfree=10G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrun.sh Cardiomyocyte

# Running 7-20-21 after upping memory
qsub -P trapnelllab -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrun.sh Fibroblast 
qsub -P trapnelllab -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrun.sh Macrophage
qsub -P trapnelllab -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrun.sh Vascular_Endothelium 
qsub -P trapnelllab -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrun.sh T_Cell


# Looks like fibroblast and macrophage were crashing. Run again with more memory...?
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrun.sh Fibroblast
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrun.sh Macrophage

qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrun.sh Macrophage

Rscript plotOutputs.R -c Cardiomyocyte


Rscript plotOutputs.R -c Endocardium
Rscript plotOutputs.R -c Vascular_Endothelium
Rscript plotOutputs.R -c T_Cell
Rscript plotOutputs.R -c Cardiomyocyte




# 7-28-21: Add age
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Fibroblast Anatomical_Site,Age
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Macrophage Anatomical_Site,Age
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Age
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Age
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Endocardium Anatomical_Site,Age

qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Log10Age
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Log10Age
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Endocardium Anatomical_Site,Log10Age
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Fibroblast Anatomical_Site,Log10Age
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Macrophage Anatomical_Site,Log10Age


qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Log10Age

# 7 -29-21 adding some plots
Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Log10Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult

Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult
Rscript plotOutputs.R -c T_Cell -m _fix_Anatomical_Site,Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult







# 7-30-21: Add age, fixed. Had a bug where I left in a part of code subsetting down to only 20 genes (meant for testing, I failed to return to original behavior of fitting all)
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Fibroblast Anatomical_Site,Age
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Macrophage Anatomical_Site,Age
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Age
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Age
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Endocardium Anatomical_Site,Age

qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Cardiomyocyte Anatomical_Site,Log10Age
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Vascular_Endothelium Anatomical_Site,Log10Age
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Endocardium Anatomical_Site,Log10Age
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Fibroblast Anatomical_Site,Log10Age
qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/singleMMrunSetFixed.sh Macrophage Anatomical_Site,Log10Age









# 8-2-21 adding some plots, fixed issue adding age
Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Log10Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult
Rscript plotOutputs.R -c Endocardium -m _fix_Anatomical_Site,Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult

Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Log10Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult
Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult

Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Log10Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult
Rscript plotOutputs.R -c Fibroblast -m _fix_Anatomical_Site,Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult

Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Log10Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult

Rscript plotOutputs.R -c T_Cell -m _fix_Anatomical_Site,Log10Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult
Rscript plotOutputs.R -c T_Cell -m _fix_Anatomical_Site,Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult




Rscript plotOutputs.R -c Vascular_Endothelium -m _fix_Anatomical_Site,Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult
Rscript plotOutputs.R -c T_Cell -m _fix_Anatomical_Site,Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult
Rscript plotOutputs.R -c Cardiomyocyte -m _fix_Anatomical_Site,Age_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult


















