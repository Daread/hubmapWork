import os
import argparse
import easygrid

# Stages
MAKE_TEMPLATE_FILES = 'make_template_files'
TRAIN_MODELS = 'train_models'
SCORE_SUMSTATS = 'score_sumstats'
GATHER_SUMSTATS = 'gather_sumstats'

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to generate LD score regression scores and plot results for a given set of genomic intervals.')
    parser.add_argument('--output_directory', default='.', help='Directory in which pipeline outputs are saved.')
    parser.add_argument('--score_sumstats_options', help='String containing options for use with LDSC in scoring sumstats files. Must start with a space and be quoted. Example: " --chisq-max 99999999999"')
    parser.add_argument('--sample_sheet', required=True, help='TSV file with headers sample_id, sites')
    parser.add_argument('--sumstats', required=True, help="TSV file with headers phenotype and sumstats for a phenotype description column and sumstats file for use with LD score regression.")
    parser.add_argument('--liftover_chain', help='Liftover chain to convert SNPs in sumstats files to the genome used for the intervals.')
    parser.add_argument('--master_peaks', help='BED file with set of peaks that make up all possible intervals that the specific set can draw from.')
    parser.add_argument('--baseline_prefix', default="/net/gs/vol1/home/ajh24/bin/ldsc/1000G_Phase3_baseline_v1.1_rerun/baseline.", help='Prefix for baseline model file to use with LD score regression.')
    parser.add_argument('--annot_template_prefix', default='/net/gs/vol1/home/ajh24/bin/ldsc/1000G_Phase3_cell_type_groups/cell_type_group.1', help='Prefix for .annot.gz files to use as templates for new annotation files.')
    parser.add_argument('--bfile_prefix', default="/net/gs/vol1/home/ajh24/bin/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC", help="Prefix for files used in --bfile argument to LD score regression.")
    parser.add_argument('--hapmap_snps_prefix', default="/net/gs/vol1/home/ajh24/bin/ldsc/hapmap3_snps/hm", help="Prefix for files used in --hapmap_snps argument to LD score regression.")
    parser.add_argument('--ld_score_prefix', default="/net/gs/vol1/home/ajh24/bin/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.", help='Prefix for --w-ld-chr argument to LD score regression.')
    parser.add_argument('--frqfile_prefix', default="/net/gs/vol1/home/ajh24/bin/ldsc/1000G_Phase3_frq/1000G.EUR.QC.", help='Prefix for --frqfile-chr argument used with LD score regression.')
    parser.add_argument('--dry', action='store_true', help='Set flag to perform dry run.')
    args = parser.parse_args()

    # Check important inputs
    easygrid.check_exists(args.sample_sheet)
    easygrid.check_exists(args.sumstats)
    easygrid.check_exists(args.master_peaks)

    if args.liftover_chain:
        easygrid.check_exists(args.liftover_chain)

    # Construct output directories
    make_template_files_dir = os.path.join(args.output_directory, MAKE_TEMPLATE_FILES)
    train_models_dir = os.path.join(args.output_directory, TRAIN_MODELS)
    score_sumstats_dir = os.path.join(args.output_directory, SCORE_SUMSTATS)
    log_dir = os.path.join(args.output_directory, '.easygrid')

    # Parse input files
    samples = list(easygrid.read_delim(args.sample_sheet))
    sumstats = list(easygrid.read_delim(args.sumstats))

    pipeline = easygrid.JobManager(log_dir)

    # Check to make sure all files specified exist
    for sample in samples:
        easygrid.check_exists(sample['sites'])

    for sumstat in sumstats:
        easygrid.check_exists(sumstat['sumstats'])

    for chromosome in range(1, 23):
        # Make a new template file for each chromosome with coordinates lifted over and such
        original_template = '%s.%s.annot.gz' % (args.annot_template_prefix, chromosome)
        easygrid.check_exists(original_template)

        template_file = os.path.join(make_template_files_dir, os.path.basename(original_template))
        command = 'python /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/ldscore_regression_scripts/generate_ld_score_annotation_template.py %s %s --peaks %s --buffer 100' % (original_template, template_file, args.master_peaks)

        if args.liftover_chain:
            command += ' --liftover_chain %s' % (args.liftover_chain)

        pipeline.add(command, name=MAKE_TEMPLATE_FILES, walltime="10:00:00", memory='10G', outputs=[template_file])

        # Now generate a model for each sample/chromosome using the templates
        for sample in samples:
            bfile = '%s.%s' % (args.bfile_prefix, chromosome)

            hapmap_file = '%s.%s.snp' % (args.hapmap_snps_prefix, chromosome)
            easygrid.check_exists(hapmap_file)

            model_prefix = os.path.join(train_models_dir, '%s.%s' % (sample['sample_id'], chromosome))
            M_file = '%s.l2.M' % (model_prefix)
            annot_file = '%s.annot.gz' % (model_prefix)
            ldscore_file = '%s.l2.ldscore.gz' % (model_prefix)
            M50_file = '%s.l2.M_5_50' % (model_prefix)

            command = 'python /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/ldscore_regression_scripts/generate_cluster_ldsc_model.py %s %s %s --bfile %s --hapmap_snps %s' % (template_file, sample['sites'], model_prefix, bfile, hapmap_file)
            pipeline.add(command, name=TRAIN_MODELS, dependencies=[MAKE_TEMPLATE_FILES], outputs=[M_file, annot_file, ldscore_file, M50_file], memory='3G', walltime='10:00:00')

    # Score each sample for every sumstats file
    results_files = []
    for sample in samples:
        model_prefix = os.path.join(train_models_dir, '%s.' % sample['sample_id'])

        for sumstat in sumstats:
            score_prefix = os.path.join(score_sumstats_dir, '%s-%s' % (sample['sample_id'], sumstat['phenotype']))
            results_file = '%s.results' % score_prefix
            results_files.append(results_file)
            # command = 'module unload python; module load python/2.7.3; python /net/gs/vol1/home/ajh24/bin/ldsc/ldsc.py --h2 %s --w-ld-chr %s --ref-ld-chr %s,%s --overlap-annot --frqfile-chr %s --out %s --print-coefficients' % (sumstat['sumstats'], args.ld_score_prefix, model_prefix, args.baseline_prefix, args.frqfile_prefix, score_prefix)
            #command = 'module unload python; module load python/2.7.13; python /net/gs/vol1/home/ajh24/bin/ldsc/ldsc.py --h2 %s --w-ld-chr %s --ref-ld-chr %s,%s --overlap-annot --frqfile-chr %s --out %s --print-coefficients' % (sumstat['sumstats'], args.ld_score_prefix, model_prefix, args.baseline_prefix, args.frqfile_prefix, score_prefix)
            command = 'module unload python; module load python/2.7.13; module load numpy/1.16.6; module load pandas/0.24.2; module load scipy/1.2.3; module load bitarray/0.8.3; module load pysam/0.16.0.1; module load pybedtools/0.7.10; python /net/gs/vol1/home/ajh24/bin/ldsc/ldsc.py --h2 %s --w-ld-chr %s --ref-ld-chr %s,%s --overlap-annot --frqfile-chr %s --out %s --print-coefficients' % (sumstat['sumstats'], args.ld_score_prefix, model_prefix, args.baseline_prefix, args.frqfile_prefix, score_prefix)
            # command = 'module unload python; module load python/2.7.13; module load numpy/1.16.6; module load pandas/0.24.2; module load scipy/1.2.3; module load bitarray/0.8.3; module load pysam/0.16.0.1; module load pybedtools/0.7.10; python /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/ldscore_regression_scripts/ldsc.py --h2 %s --w-ld-chr %s --ref-ld-chr %s,%s --overlap-annot --frqfile-chr %s --out %s --print-coefficients' % (sumstat['sumstats'], args.ld_score_prefix, model_prefix, args.baseline_prefix, args.frqfile_prefix, score_prefix)
 

            # Add DFR 9-8-21
            print(command)

            if args.score_sumstats_options:
                command += args.score_sumstats_options
            pipeline.add(command, name=SCORE_SUMSTATS, dependencies=[TRAIN_MODELS], memory='8G', outputs=[results_file], walltime="10:00:00")

    # Gather the results file stats
    results_file_list = os.path.join(args.output_directory, 'results_file_list.txt')
    results_file_gathered = os.path.join(args.output_directory, 'results_gathered.txt')

    command = 'python /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/ldscore_regression_scripts/gather_score_sumstats_results.py %s %s' % (results_file_list, results_file_gathered)
    pipeline.add(command, name=GATHER_SUMSTATS, dependencies=[SCORE_SUMSTATS], memory='5G', outputs=[results_file_gathered], walltime='10:00:00')

    # Everything seems ok, so set up directories and run
    if not args.dry:
        easygrid.mkdir(args.output_directory)
        easygrid.mkdir(make_template_files_dir)
        easygrid.mkdir(train_models_dir)
        easygrid.mkdir(score_sumstats_dir)

    with open(results_file_list, 'w') as results_out:
        for result_file in results_files:
            results_out.write(result_file + '\n')

    pipeline.run(infer_dependencies=False, dry=args.dry)
