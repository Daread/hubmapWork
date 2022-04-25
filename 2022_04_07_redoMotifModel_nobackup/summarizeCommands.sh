

# Looking at promoter-only models:


Rscript ./summarizeFitResults.R --featureSelection Binary_PromOnly_Motif_Counts --summarizationPolicy alpha0.5_max --predictionFraming regression --paramFile VaryPromoters
Rscript ./summarizeFitResults.R --featureSelection Binary_PromOnly_Motif_Counts --summarizationPolicy alpha0.5_average --predictionFraming regression --paramFile VaryPromoters

Rscript ./summarizeFitResults.R --featureSelection Binary_PromOnly_Motif_Counts --summarizationPolicy alpha0.5_max --predictionFraming classification --paramFile VaryPromoters
Rscript ./summarizeFitResults.R --featureSelection Binary_PromOnly_Motif_Counts --summarizationPolicy alpha0.5_average --predictionFraming classification --paramFile VaryPromoters






# 9-22-21 Final fitting sets of runs
Rscript ./summarizeFitResults.R --featureSelection Binary_Combined_Motif_Counts --summarizationPolicy alpha0.5_max --predictionFraming regression --paramFile VaryDistalParams
Rscript ./summarizeFitResults.R --featureSelection Binary_Combined_Motif_Counts --summarizationPolicy alpha0.5_average --predictionFraming regression --paramFile VaryDistalParams

Rscript ./summarizeFitResults.R --featureSelection Binary_Combined_Motif_Counts --summarizationPolicy alpha0.5_max --predictionFraming classification --paramFile VaryDistalParams
Rscript ./summarizeFitResults.R --featureSelection Binary_Combined_Motif_Counts --summarizationPolicy alpha0.5_average --predictionFraming classification --paramFile VaryDistalParams






# 9-22-21 Final fitting sets of runs
Rscript ./summarizeFitResults.R --featureSelection Binary_Combined_Motif_Counts --summarizationPolicy alpha0.5_max --predictionFraming regression --paramFile VaryDistalAllRunsMade
Rscript ./summarizeFitResults.R --featureSelection Binary_Combined_Motif_Counts --summarizationPolicy alpha0.5_average --predictionFraming regression --paramFile VaryDistalAllRunsMade
Rscript ./summarizeFitResults.R --featureSelection Binary_Combined_Motif_Counts --summarizationPolicy alpha0.5_max --predictionFraming classification --paramFile VaryDistalAllRunsMade
Rscript ./summarizeFitResults.R --featureSelection Binary_Combined_Motif_Counts --summarizationPolicy alpha0.5_average --predictionFraming classification --paramFile VaryDistalAllRunsMade




