#!/usr/bin/python
import sys
import os
import argparse
import pandas as pd
import re
import random
import string


global VCF_COLS
VCF_COLS = {}
VCF_COLS['VariantCalling']= ["CHROM","POS","ID","REF","ALT","QUAL",
							"FILTER","INFO","FORMAT","SAMPLE"]
VCF_COLS['VariantAnnotation']= ["Uploaded_variation","Location","Allele",
								"Gene","Feature","Feature_type",
								"Consequence","cDNA_position","CDS_position",
								"Protein_position","Amino_acids","Codons",
								"Existing_variation","Extra"]


def parse_args(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-t', '--tool',
						help='Tool for variation calling or variant annotation. Choose from=%s' % [k for k in VCF_COLS.keys()])
	parser.add_argument('-i', '--input_vcf', 
						help='Filepath of vcf as input.')
	parser.add_argument('-o', '--output_vcf', 
						help='Filepath of the filtered vcf.')
	parser.add_argument('-g', '--gene_regions',
						help='Filepath of the gene-of-interest regions in bed format.')
	return parser.parse_args(argv[1:])


def load_vcf(filepath, tool):
	## Load vcf and split into header dataframe and vcf dataframe
	vcf = pd.read_csv(filepath, names=VCF_COLS[tool], sep='\t')
	for i,x in vcf[VCF_COLS[tool][0]].iteritems():
		if x.startswith('#'+VCF_COLS[tool][0]):
			header_indx = i
			break
	header = vcf.loc[header_indx,].to_frame().T
	vcf = vcf.loc[range(header_indx+1, vcf.shape[0]),]
	# vcf = vcf.reset_index()
	return (header, vcf)


def filter_multi_alt(vcf):
	## Remove site with >1 alternative allele
	return vcf.loc[[',' not in x for x in vcf['ALT']], ]


def filter_mutect_info_format(vcf):
	## Filter INFO and FORMAT field based on the threshold
	filter_name = 'DPl30-ECNTl5-AFg95-l55-g45-l05'
	for i,row in vcf.iterrows():
		attrs = {}
		attrs['DP'] = get_info_attribute(row['INFO'], 'DP')
		attrs['ECNT'] = get_info_attribute(row['INFO'], 'ECNT')
		attrs['AF'] = get_format_attribute(row['FORMAT'], row['SAMPLE'], 'AF')
		## apply thresholding
		if attrs['DP'] <= 30:
			vcf.loc[i, 'FILTER'] = filter_name
		elif attrs['ECNT'] <=5:
			vcf.loc[i, 'FILTER'] = filter_name
		elif attrs['AF'] >= .95 or attrs['AF'] <= 0.05 or (attrs['AF'] <= 0.55 and attrs['AF'] >= 0.45):
			vcf.loc[i, 'FILTER'] = filter_name
		else:
			vcf.loc[i, 'FILTER'] = 'PASS'
	return vcf


def filter_varscan_info_format(vcf):
	## Filter INFO and FORMAT field based on the threshold
	filter_name = 'ADl6-SDPl20-ADFle0-ADRle0-FREQl05'
	for i,row in vcf.iterrows():
		attrs = {}
		attrs['AD'] = get_format_attribute(row['FORMAT'], row['SAMPLE'], 'AD')
		attrs['SDP'] = get_format_attribute(row['FORMAT'], row['SAMPLE'], 'SDP')
		attrs['FREQ'] = get_format_attribute(row['FORMAT'], row['SAMPLE'], 'FREQ')
		attrs['ADF'] = get_format_attribute(row['FORMAT'], row['SAMPLE'], 'ADF')
		attrs['ADR'] = get_format_attribute(row['FORMAT'], row['SAMPLE'], 'ADR')
		## apply thresholding
		# if (attrs['AD'] <= 6) or (attrs['SDP'] <= 20) or (attrs['FREQ'] >= .95 or attrs['FREQ'] <= 0.05 or (attrs['FREQ'] <= 0.55 and attrs['FREQ'] >= 0.45)):
		if (attrs['AD'] <= 6) or (attrs['SDP'] <= 20) or attrs['FREQ'] >= 0.45 or (attrs['ADF'] <= 0 or attrs['ADR'] <= 0):
			vcf.loc[i, 'FILTER'] = filter_name
		else:
			vcf.loc[i, 'FILTER'] = 'PASS'
	return vcf


def get_info_attribute(field, attr):
	matches = re.findall(r"%s=(.+?);" % attr, field)
	if len(matches) > 1:
		sys.exit("ERROR: Should not have duplicate attribute %s in FORMAT." % attr)
	elif len(matches) < 1:
		sys.exit("ERROR: No attribute %s found in INFO." % attr)
	else:
		return float(matches[0])


def get_format_attribute(header, field, attr):
	header = header.split(':')
	indx = [i for i in range(len(header)) if header[i]==attr]
	if len(indx) > 1:
		sys.exit("ERROR: Should not have duplicate attribute %s in FORMAT." % attr)
	elif len(indx) < 1:
		sys.exit("ERROR: No attribute %s found in FORMAT." % attr)
	else:
		value = field.split(':')[indx[0]]
		value = float(value.strip('%'))/100 if value.endswith('%') else float(value)
		return value


def filter_gene_region(filepath_in, region):
	## Narrow down regions to the genes of interest
	filepath_tmp = '/tmp/'+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))+'.vcf'
	os.system('grep "#" '+ filepath_in +' > '+ filepath_tmp)
	os.system('bedtools intersect -a '+filepath_in +' -b '+ region +' >> '+ filepath_tmp)
	header, vcf = load_vcf(filepath_tmp, "VariantCalling")
	os.system('rm '+ filepath_tmp)
	return vcf


def filter_vep_tumor_artifact(vcf):
	## Filter tumor artifacts called by VEP 
	filter = ["upstream_gene_variant",
			"downstream_gene_variant",
			"non_coding_transcript_exon_variant",
			"non_coding_transcript_variant",
			"synonymous_variant",
			"intron_variant"]
	valid_indx = []
	for i,row in vcf.iterrows():
		if all([x not in filter for x in row["Consequence"].split(',')]):
			valid_indx.append(i)
	return vcf.loc[valid_indx,:]


def save_filtered_vcf(filepath_in, filepath_out, vcf):
	## Save vcf file by adding back the header
	filepath_tmp = '/tmp/'+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))+'.vcf'
	vcf.to_csv(filepath_tmp, sep='\t', header=False, index=False)
	os.system('grep "#" '+ filepath_in +' > '+ filepath_out)
	os.system('cat '+ filepath_tmp +' >> '+ filepath_out)
	os.system('rm '+ filepath_tmp)


def main(argv):
	parsed = parse_args(argv)
	## filter variant calling output
	if parsed.tool == "VariantCalling":
		tmp_vcf = parsed.output_vcf+'.tmp'
		_, vcf = load_vcf(parsed.input_vcf, parsed.tool)
		vcf = filter_multi_alt(vcf)
		vcf = filter_varscan_info_format(vcf)
		save_filtered_vcf(parsed.input_vcf, tmp_vcf, vcf)
		vcf = filter_gene_region(tmp_vcf, parsed.gene_regions)
		save_filtered_vcf(parsed.input_vcf, parsed.output_vcf, vcf)
		os.system('rm '+ tmp_vcf)
	## filter var 
	elif parsed.tool == "VariantAnnotation":
		_, vcf = load_vcf(parsed.input_vcf, parsed.tool)
		vcf = filter_vep_tumor_artifact(vcf)
		save_filtered_vcf(parsed.input_vcf, parsed.output_vcf, vcf)
	else:
		sys.exit("ERROR: Choose from %s" % [k for k in VCF_COLS.keys()])

if __name__ == '__main__':
	main(sys.argv)
