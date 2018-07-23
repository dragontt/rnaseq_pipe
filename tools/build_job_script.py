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
	parser.add_argument('-f', '--reference_genome', 
						help='Annotation reference file in GTF/GFF3 format.')
	# parser.add_argument('-g', '--gene_list',
	# 					help='Gene list.')
	return parser.parse_args(argv[1:])


def build_alignment(samples, genome_index):
	built_dict = {}
	for i,row in samples.iterrows():
		## define variables
		sample = '_'.join([row['RUN'],row['SAMPLE']])
		seqs = [x.strip('" ') for x in row['SEQUENCE'].split(',')]
		paired = True if len(seqs) == 2 else False
		target_dir = 'alignment/'+ sample
		bam_file = target_dir +'/tmp.bam'
		sorted_bam_file = target_dir +'/alignment.bam'
		aligned_bam_file = target_dir +'/aligned_only.bam'
		log_file = target_dir +'/alignment_summary.log'
		## create cleaned directory
		recipe = 'rm -rf '+ target_dir +'; mkdir '+ target_dir +'; '
		## align reads
		recipe += 'hisat2 -p 10 --rg-id='+ sample +'--rg SM:'+ sample +' -x '+ genome_index
		recipe += ' -1 '+ seqs[0] +' -2 '+ seqs[1] if paried else ' -U '+ seqs[0]
		recipe += ' -S '+ bam_file +' 2> '+ log_file +'; '
		## sort bam file and index
		recipe += 'samtools sort -@ 10 -T /home/tmp/'+ sample +'.sort_bam -o '+ sorted_bam_file +' '+ bam_file +'; '
		recipe += 'rm ' + bam_file
		recipe += 'samtools index'+ sorted_bam_file +'; '
		## get bam of only aligned reads
		recipe += 'samtools view -b -f 2 '+ sorted_bam_file +' > '+ aligned_bam_file +'; '
		recipe += 'samtools index'+ aligned_bam_file +'; '
		## set alignment dictionary
		built_dict[aligned_bam_file] = {'prereq':' '.join(seqs), 'recipe':recipe}
	return built_dict


def build_expression_quantification(samples, ref_gtf, tool='htseq'):
	if tool not in ['htseq', 'stringtie']:
		sys.exit('ERROR: Wrong expression quantification tool specified.')
	built_dict = {}
	for i,row in samples.iterrows():
		sample = '_'.join([row['RUN'],row['SAMPLE']])
		prereq = 'alignment/hisat2/'+ sample +'/alignment.bam'
		## use htseq-count for raw count
		if tool == 'htseq':
			target_dir = 'expression/htseq/'+ sample
			output_file = target_dir +'/gene_count.tsv'
			recipe = 'rm -rf '+ target_dir +'; mkdir '+ target_dir +'; '
			recipe += 'htseq-count -f bam -i gene_id -s yes -t gene '+ prereq +' '+ ref_gtf +' > '+ output_file
		## use stringtie for fpkm/tpm
		elif tool == 'stringtie':
			target_dir = 'expression/stringtie/'+ sample
			output_gtf = target_dir +'/gene_expression.gtf'
			output_file = target_dir +'/gene_abundances.tab'
			recipe = 'rm -rf '+ target_dir +'; mkdir '+ target_dir +'; '
			recipe += 'stringtie '+ prereq +' -G '+ ref_gtf +' -e -o '+ output_gtf +' -A '+ output_file
		## set expression dict
		built_dict[output_file] = {'prereq': prereq, 'recipe': recipe}
	return built_dict


def build_variant_calling(samples, ref_genome, xxx, tool="mutect2"):
	built_dict = {}
	for i,row in samples.iterrows():
		sample = '_'.join([row['RUN'],row['SAMPLE']])
		prereq = 'alignment/hisat2/'+ sample +'/aligned_only.bam'

		target_dir = 'var_calling/'+ tool +'/'+ sample
		output_vcf = target_dir +'/tumor_variants.vcf.gz'
		log_prefix = target_dir +'/'+ tool
		recipe = 'rm -rf '+ target_dir +'; mkdir '+ target_dir +'; '
		recipe += 'gatk Mutect2 -R '+ ref_genome +' -I '+ prereq +' -tumor '+ sample +' -VS STRICT --germline-resource '+ xxx +'--af-of-alleles-not-in-resource 0.0000025 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -O '+ output_vcf +' > '+ log_prefix +'.out 2> '+ log_prefix +'.err;'
		## set vc dict
		built_dict[output_file] = {'prereq': prereq, 'recipe': recipe}
	return built_dict


def main(argv):
	parsed = parse_args(argv)
	samples = pd.read_csv(parsed.samples, delimiter='\t', dtype=str)
	
	align_dict = build_alignment(samples, parsed.genome_index)
	for target in sorted(align_dict.keys()):
		sys.stdout.write('%s: %s\n\t%s\n' % (target, align_dict[target]['prereq'], align_dict[target]['recipe']))
	sys.stdout.write('ALIGNMENTS = %s\n' % ' '.join(sorted(align_dict.keys())))

	# sys.stdout.write('all: $(ALIGNMENTS)')

	expr_dict = build_expression_quantification(samples, parsed.reference_gtf)
	for target in sorted(expr_dict.keys()):
		sys.stdout.write('%s: %s\n\t%s\n' % (target, expr_dict[target]['prereq'], expr_dict[target]['recipe']))
	sys.stdout.write('EXPRESSIONS = %s\n' % ' '.join(sorted(expr_dict.keys())))

	# sys.stdout.write('all: $(EXPRESSIONS)')

	vc_dict = build_variant_calling(samples, parsed.reference_genome, xxx)
	for target in sorted(vc_dict.keys()):
		sys.stdout.write('%s: %s\n\t%s\n' % (target, vc_dict[target]['prereq'], vc_dict[target]['recipe']))
	sys.stdout.write('VARCALLS = %s\n' % ' '.join(sorted(expr_dict.keys())))

	# sys.stdout.write('all: $(VARCALLS)')

	sys.stdout.write('all: $(ALIGNMENTS) $(EXPRESSIONS) $(VARCALLS)')

if __name__ == '__main__':
	main(sys.argv)
