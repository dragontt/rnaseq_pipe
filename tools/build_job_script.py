#!/usr/bin/python
import sys
import argparse
import pandas as pd


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--samples', 
						help='Sample summary metadata file.')
	parser.add_argument('-i', '--genome_index', 
						help='Prefix of genome index files built for the aligner, e.g. hisat2')
	parser.add_argument('-r', '--reference_gtf', 
						help='Annotation reference file in GTF/GFF3 format.')
	parser.add_argument('-g', '--gene_list',
						help='Gene list.')
	return parser.parse_args(argv[1:])


def build_alignment(samples, genome_index):
	## define alignment dictionary
	align_dict = {}
	for i,row in samples.iterrows():
		## define variables
		sample = row['SAMPLE']
		seqs = [x.strip('" ') for x in row['SEQUENCE'].split(',')]
		paired = True if len(seqs) == 2 else False
		target_dir = 'alignment/hisat2/'+ sample
		bam_file = target_dir +'/aligned_reads.bam'
		sorted_bam_file = target_dir +'/aligned_reads_sorted.bam'
		log_file = target_dir +'/alignment_summary.log'
		## create cleaned directory
		recipe = 'rm -rf '+ target_dir +'; mkdir '+ target_dir +'; '
		## align reads
		## TODO: add arugment for trim5
		if paired:
			recipe += 'hisat2 -p 2 --rg-id='+ sample + ' -x '+ genome_index +' --dta --rna-strandness RF -1 '+ seqs[0] +' -2 '+ seqs[1] +' -S '+ bam_file +' 2> '+ log_file +'; '
		else:
			recipe += 'hisat2 -p 2 --rg-id='+ sample + ' -x '+ genome_index +' --dta --rna-strandness RF -U '+ seqs[0] +' -S '+ bam_file +' 2> '+ log_file +'; '
		## sort bam file
		recipe += 'samtools sort -T /tmp/'+ sample +'.sort_bam -o '+ sorted_bam_file +' '+ bam_file +'; '
		recipe += 'rm ' + bam_file
		## set alignment dictionary
		align_dict[sorted_bam_file] = {'prereq': ' '.join(seqs), 
										'recipe': recipe}
	return align_dict


def build_expression_quantification(samples, reference_gtf):
	expr_dict = {}
	for i,row in samples.iterrows():
		sample = row['SAMPLE']
		prereq = 'alignment/hisat2/'+ sample +'/aligned_reads_sorted.bam'
		target_dir = 'expression/stringtie/'+ sample
		output_gtf_file = target_dir +'/gene_expression.gtf'
		output_tab_file = target_dir +'/gene_abundances.tab'
		recipe = 'rm -rf '+ target_dir +'; mkdir '+ target_dir +'; '
		recipe += 'stringtie '+ prereq +' -G '+ reference_gtf +' -e -o '+ output_gtf_file +' -A '+ output_tab_file
		expr_dict[output_tab_file] = {'prereq': prereq, 
										'recipe': recipe}
	return expr_dict


def main(argv):
	parsed = parse_args(argv)
	samples = pd.read_csv(parsed.samples, delimiter='\t', dtype=str)
	
	align_dict = build_alignment(samples, parsed.genome_index)
	for target in sorted(align_dict.keys()):
		sys.stdout.write('%s: %s\n\t%s\n' % (target, align_dict[target]['prereq'], align_dict[target]['recipe']))
	sys.stdout.write('HISAT2_ALIGNMENTS = %s\n' % ' '.join(sorted(align_dict.keys())))

	expr_dict = build_expression_quantification(samples, parsed.reference_gtf)
	for target in sorted(expr_dict.keys()):
		sys.stdout.write('%s: %s\n\t%s\n' % (target, expr_dict[target]['prereq'], expr_dict[target]['recipe']))
	sys.stdout.write('STRINGTIE_EXPRESSIONS = %s\n' % ' '.join(sorted(expr_dict.keys())))

	sys.stdout.write('all: $(HISAT2_ALIGNMENTS) $(STRINGTIE_EXPRESSIONS)')

if __name__ == '__main__':
	main(sys.argv)
