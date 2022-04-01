# Data from an Asthma GWAS called gabriel
# Note https://beaune.cng.fr/gabriel/gabriel_results_description.xls is a description of the columns
wget https://beaune.cng.fr/gabriel/gabriel_results.zip
unzip gabriel_results.zip
rm gabriel_results.zip
mv gabriel_asthma_meta-analysis_36studies_format_repository_NEJM.txt gabriel.asthma.txt 

## Chronic Kidney disease downloaded from NHLBI (http://fox.nhlbi.nih.gov/)
wget http://fox.nhlbi.nih.gov/CKDGen/formatted_round3meta_CKD_overall_IV_2GC_b36_MAFget005_Nget50_20120725_b37.csv.gz
zcat formatted_round3meta_CKD_overall_IV_2GC_b36_MAFget005_Nget50_20120725_b37.csv.gz | sed 's/,/\t/g' > ckdgen.chronic_kidney_disease.txt
rm formatted_round3meta_CKD_overall_IV_2GC_b36_MAFget005_Nget50_20120725_b37.csv.gz

## IgA defficiency Bronson et al Nat. Gen from EMBL summary stats
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/BronsonPG_27723758/BronsonEtAl_NatGenet_2016.zip
unzip BronsonEtAl_NatGenet_2016.zip

mv BronsonEtAl_NatGenet_2016/Bronson_et_al.2016.IgAD_GWAS_Meta_Analysis_Summary_Stats.txt bronson.iga_deficiency.txt
rm -r BronsonEtAl_NatGenet_2016*

# Reproductive stuff from Barban Nature Genetics from EMBL summary stats
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/BarbanN_27798627/NumberChildrenEverBorn_Male.txt.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/BarbanN_27798627/AgeFirstBirth_Male.txt.gz

zcat NumberChildrenEverBorn_Male.txt.gz > barban.number_children_born_male.txt
zcat AgeFirstBirth_Male.txt.gz > barban.age_first_birth_male.txt

rm NumberChildrenEverBorn_Male.txt.gz
rm AgeFirstBirth_Male.txt.gz

# These studies were all downloaded from https://data.broadinstitute.org/alkesgroup/sumstats_formatted/ (linked to from Broad LD Hub)
# They are already sumstats formatted so don't need to munge them
## Multiple Sclerosis from https://www.ncbi.nlm.nih.gov/pubmed/21833088
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Multiple_sclerosis.sumstats
cat PASS_Multiple_sclerosis.sumstats | gzip > ldsc_sumstats_files/imsgc.multiple_sclerosis.sumstats.gz
rm PASS_Multiple_sclerosis.sumstats

## Rheumatoid Arthritis from https://www.nature.com/nature/journal/v506/n7488/full/nature12873.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Rheumatoid_Arthritis.sumstats
cat PASS_Rheumatoid_Arthritis.sumstats | gzip > ldsc_sumstats_files/okada.rheumatoid_arthritis.sumstats.gz
rm PASS_Rheumatoid_Arthritis.sumstats

## Type 1 diabetes from http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002293
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Type_1_Diabetes.sumstats
cat PASS_Type_1_Diabetes.sumstats | gzip > ldsc_sumstats_files/bradfield.type1_diabetes.sumstats.gz
rm PASS_Type_1_Diabetes.sumstats

## Type 2 diabetes from http://www.nature.com/ng/journal/v44/n9/full/ng.2383.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Type_2_Diabetes.sumstats
cat PASS_Type_2_Diabetes.sumstats | gzip > ldsc_sumstats_files/morris.type2_diabetes.sumstats.gz
rm PASS_Type_2_Diabetes.sumstats

## Celiac Disease from https://www.ncbi.nlm.nih.gov/pubmed/20190752
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Celiac.sumstats
cat PASS_Celiac.sumstats | gzip > ldsc_sumstats_files/dubois.celiac_disease.sumstats.gz
rm PASS_Celiac.sumstats

## Lupus from http://www.nature.com/ng/journal/v47/n12/abs/ng.3434.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Lupus.sumstats
cat PASS_Lupus.sumstats | gzip >  ldsc_sumstats_files/bentham.lupus.sumstats.gz
rm PASS_Lupus.sumstats

## Height from https://www.nature.com/nature/journal/v467/n7317/abs/nature09410.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Height1.sumstats
cat PASS_Height1.sumstats | gzip > ldsc_sumstats_files/lango_allen.height.sumstats.gz
rm PASS_Height1.sumstats 

## Fasting glucose levels from https://www.ncbi.nlm.nih.gov/pubmed/22581228
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Fasting_Glucose.sumstats
cat PASS_Fasting_Glucose.sumstats | gzip > ldsc_sumstats_files/manning.fasting_glucose.sumstats.gz
rm PASS_Fasting_Glucose.sumstats

## BMI from https://www.nature.com/ng/journal/v42/n11/abs/ng.686.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_BMI1.sumstats
cat PASS_BMI1.sumstats | gzip > ldsc_sumstats_files/speliotes.bmi.sumstats.gz
rm PASS_BMI1.sumstats

## Smoking behavior from http://www.nature.com/ng/journal/v42/n5/abs/ng.571.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Ever_Smoked.sumstats
cat PASS_Ever_Smoked.sumstats | gzip > ldsc_sumstats_files/tgc.ever_smoked.sumstats.gz
rm PASS_Ever_Smoked.sumstats

## Alzheimers disease http://www.nature.com/ng/journal/v45/n12/full/ng.2802.html 
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Alzheimer.sumstats
cat PASS_Alzheimer.sumstats | gzip > ldsc_sumstats_files/lambert.alzheimers.sumstats.gz
rm PASS_Alzheimer.sumstats

## Schizophrenia from https://www.nature.com/nature/journal/v511/n7510/full/nature13595.html 
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Schizophrenia.sumstats
cat PASS_Schizophrenia.sumstats | gzip > ldsc_sumstats_files/pgc.schizophrenia.sumstats.gz
rm PASS_Schizophrenia.sumstats

## Bipolar disorder from https://www.ncbi.nlm.nih.gov/pubmed/21926972
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Bipolar_Disorder.sumstats
cat PASS_Bipolar_Disorder.sumstats | gzip > ldsc_sumstats_files/pgc.bipolar_disorder.sumstats.gz
rm PASS_Bipolar_Disorder.sumstats

## Autism from ross disorder study http://www.thelancet.com/journals/lancet/article/PIIS0140-6736(12)62129-1/abstract
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Autism.sumstats
cat PASS_Autism.sumstats | gzip > ldsc_sumstats_files/pgc.autism.sumstats.gz
rm PASS_Autism.sumstats

## Educational attainment from http://science.sciencemag.org/content/340/6139/1467
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Years_of_Education1.sumstats
cat PASS_Years_of_Education1.sumstats | gzip > ldsc_sumstats_files/rietveld.educational_attainment.sumstats.gz
rm PASS_Years_of_Education1.sumstats

## More educational attainment from http://www.nature.com/nature/journal/v533/n7604/full/nature17671.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Years_of_Education2.sumstats
cat PASS_Years_of_Education2.sumstats | gzip > ldsc_sumstats_files/okbay.educational_attainment.sumstats.gz
rm PASS_Years_of_Education2.sumstats

## Neuroticism from https://www.nature.com/ng/journal/v48/n6/full/ng.3552.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Neuroticism.sumstats
cat PASS_Neuroticism.sumstats | gzip > ldsc_sumstats_files/okbay.neuroticism.sumstats.gz
rm PASS_Neuroticism.sumstats

## Self-reported well-being also from https://www.nature.com/ng/journal/v48/n6/full/ng.3552.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_SWB.sumstats
cat PASS_SWB.sumstats | gzip > ldsc_sumstats_files/okbay.self_reported_well_being.sumstats.gz
rm PASS_SWB.sumstats

## Depressive symptoms also from https://www.nature.com/ng/journal/v48/n6/full/ng.3552.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_DS.sumstats
cat PASS_DS.sumstats | gzip > ldsc_sumstats_files/okbay.depressive_symptoms.sumstats.gz
rm PASS_DS.sumstats

## HDL from https://www.nature.com/nature/journal/v466/n7307/full/nature09270.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_HDL.sumstats
cat PASS_HDL.sumstats | gzip > ldsc_sumstats_files/teslovich.hdl_cholesterol.sumstats.gz
rm PASS_HDL.sumstats

## LDL also from https://www.nature.com/nature/journal/v466/n7307/full/nature09270.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_LDL.sumstats
cat PASS_LDL.sumstats | gzip > ldsc_sumstats_files/teslovich.ldl_cholesterol.sumstats.gz
rm PASS_LDL.sumstats

## Triglycerides also from https://www.nature.com/nature/journal/v466/n7307/full/nature09270.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Triglycerides.sumstats
cat PASS_Triglycerides.sumstats | gzip > ldsc_sumstats_files/teslovich.triglycerides.sumstats.gz
rm PASS_Triglycerides.sumstats

## IBD from https://www.nature.com/nature/journal/v491/n7422/full/nature11582.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_IBD.sumstats
cat PASS_IBD.sumstats | gzip > ldsc_sumstats_files/jostins.ibd.sumstats.gz
rm PASS_IBD.sumstats

## Crohns also from https://www.nature.com/nature/journal/v491/n7422/full/nature11582.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Crohns_Disease.sumstats
cat PASS_Crohns_Disease.sumstats | gzip > ldsc_sumstats_files/jostins.crohns.sumstats.gz
rm PASS_Crohns_Disease.sumstats

## Ulcerative colitis also from https://www.nature.com/nature/journal/v491/n7422/full/nature11582.html
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Ulcerative_Colitis.sumstats
cat PASS_Ulcerative_Colitis.sumstats | gzip > ldsc_sumstats_files/jostins.ulcerative_colitis.sumstats.gz
rm PASS_Ulcerative_Colitis.sumstats

## Primary biliary cirrhosis from  https://www.nature.com/articles/ncomms9019
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Primary_biliary_cirrhosis.sumstats
cat PASS_Primary_biliary_cirrhosis.sumstats | gzip > ldsc_sumstats_files/cordell.pb_cirrhosis.sumstats.gz
rm PASS_Primary_biliary_cirrhosis.sumstats

## Coronary artery disease from https://www.ncbi.nlm.nih.gov/pubmed/21378990
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Coronary_Artery_Disease.sumstats
cat PASS_Coronary_Artery_Disease.sumstats | gzip > ldsc_sumstats_files/schunkert.coronary_artery_disease.sumstats.gz
rm PASS_Coronary_Artery_Disease.sumstats

## Anorexia from https://www.ncbi.nlm.nih.gov/pubmed/24514567
wget https://data.broadinstitute.org/alkesgroup/sumstats_formatted/PASS_Anorexia.sumstats
cat PASS_Anorexia.sumstats | gzip > ldsc_sumstats_files/boraska.anorexia.sumstats.gz
rm PASS_Anorexia.sumstats

# Finally for any files that were not already in sumstats format, munge them with LDSC
bash run_prepare_sumstats_files.sh

# Make a metadata table
python make_sumstats_metadata.py
