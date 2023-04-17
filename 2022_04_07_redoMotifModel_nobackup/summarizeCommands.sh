

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








# 10-24-22 Look at Lyon for setting up coaccessibility links
Rscript ./summarizeFitResults.R --featureSelection Binary_Combined_Motif_Counts --summarizationPolicy alpha0.5_average --predictionFraming regression --paramFile VaryDistalLyon2







# 2-3-23 Double-checking on generation of promoter-only plot



Rscript ./submissionVersionSummarizeFitResults.R --featureSelection Binary_PromOnly_Motif_Counts --predictionFraming regression --paramFile VaryPromotersAllRunsMade






