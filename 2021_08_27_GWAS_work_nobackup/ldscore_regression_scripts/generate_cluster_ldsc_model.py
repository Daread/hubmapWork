
print("Loading modules")

import argparse
import gzip
import pandas as pd
import math
import subprocess


def parse_specific_peaks(file_object):
    columns = ['chr', 'start', 'stop']
    specific_peaks = []

    for line in file_object:
        entries = line.strip().split('\t')
        entry_dict = dict(zip(columns, entries))
        specific_peaks.append(entry_dict)

    final_dict = {}
    for peak in specific_peaks:
        final_dict['%s.%s.%s' % (peak['chr'], peak['start'], peak['stop'])] = peak

    return final_dict


def parse_template_file(file_object):
    results = []

    file_object.readline()
    for line in file_object:
        columns = ['CHR', 'BP', 'SNP', 'CM', 'LIFTOVER', 'INTERSECTS_PEAK', 'PEAK_CHR', 'PEAK_START', 'PEAK_STOP']

        entries = line.strip().split('\t')
        entries_dict = dict(zip(columns, entries))

        results.append(entries_dict)
    return results


def add_categories(template_snps, specific_peaks):
    new_peaks = []

    for snp in template_snps:

        if snp['INTERSECTS_PEAK'] == "1":
            peak = '%s.%s.%s' % (snp['PEAK_CHR'], snp['PEAK_START'], snp['PEAK_STOP'])

            if peak in specific_peaks:
                snp['CATEGORY'] = 0

        new_peaks.append(snp)

    return new_peaks


def ldsc(ldsc_path, bfile, annotation_file, output_file_prefix, hapmap_snps):
    # command = "module unload python; module load python/2.7.3; python %s --l2 --bfile %s --ld-wind-cm 1 --annot %s --out %s --print-snps %s" % (ldsc_path, bfile, annotation_file, output_file_prefix, hapmap_snps)
    command = "module unload python; module load python/2.7.13; module load numpy/1.16.6; module load pandas/0.24.2; module load scipy/1.2.3; module load bitarray/0.8.3; module load pysam/0.16.0.1; module load pybedtools/0.7.10; python %s --l2 --bfile %s --ld-wind-cm 1 --annot %s --out %s --print-snps %s" % (ldsc_path, bfile, annotation_file, output_file_prefix, hapmap_snps)

    
    # command = "module unload python; module load python/3.7.7; python %s --l2 --bfile %s --ld-wind-cm 1 --annot %s --out %s --print-snps %s" % (ldsc_path, bfile, annotation_file, output_file_prefix, hapmap_snps)
    print(command)

    subprocess.call(command, shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to take a set of cluster specific peaks and a LDSC template and generate final annotation file and LDSC model.')
    parser.add_argument('template_file', help='LDSC template file generated by generate_ld_score_annotation_template.py')
    parser.add_argument('specific_peaks_file', help='BED file with chrom, start, and end (no header row). Any other fields are ignored.')
    parser.add_argument('output_file_prefix', help='Prefix to use for LDSC output files.')
    parser.add_argument('--ldsc_path', default="/net/gs/vol1/home/ajh24/bin/ldsc/ldsc.py", help='Path to LDSC python script.')
    parser.add_argument('--bfile', help='Prefix to PLINK files from LDSC')
    parser.add_argument('--hapmap_snps', help='hapmap SNPs file to pass to LDSC')
    args = parser.parse_args()

    # DFR Add
    print("Starting generate_cluster_ldsc_model")

    # Load peaks
    specific_peaks = parse_specific_peaks(open(args.specific_peaks_file))

    # Load template SNPs and annotate with categories
    template_snps = parse_template_file(gzip.open(args.template_file, 'rt'))
    template_snps_with_categories = add_categories(template_snps, specific_peaks)

    # Output file
    output_file = '%s.annot.gz' % args.output_file_prefix
    output_columns = ['CHR', 'BP', 'SNP', 'CM']
    category_names = ['SPECIFIC_PEAK_INTERSECT']

    with gzip.open(output_file, 'wt') as out:
        out.write('\t'.join(output_columns + category_names) + '\n')

        for snp in template_snps_with_categories:
            output_entries = [snp[column] for column in output_columns]

            # Get category indicators
            category_indicator = [0]

            if 'CATEGORY' in snp:
                category_indicator[snp['CATEGORY']] = 1

            final_entries = [str(x) for x in output_entries + category_indicator]
            out.write('\t'.join(final_entries) + '\n')

    # Run LDSC on prepared file to generate model
    ldsc(args.ldsc_path, args.bfile, output_file, args.output_file_prefix, args.hapmap_snps)
