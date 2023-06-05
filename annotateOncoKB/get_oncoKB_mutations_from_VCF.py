#!/usr/bin/env python3

### extract nonsynonymous mutations in OncoKB transcripts from VCF for OncoKB annotation ###
### Ben Kinnersley (ben.kinnersley@icr.ac.uk)
### Command line: python get_oncoKB_mutations_from_VCF.py <input_vcf> <subtype> <sample_id> <transcripts_file> <token> <output_dir> <output>

import os
import sys
import gzip
import re
import requests
import json
import time
import csv

if len(sys.argv) == 8:
	input_vcf = sys.argv[1]
	subtype = sys.argv[2]
	plate_tumour = sys.argv[3]
	transcripts_file = sys.argv[4]
	token = sys.argv[5]
	output_dir = sys.argv[6]
	output = sys.argv[7]
else:
	print('insufficient arguments provided - python get_oncoKB_mutations_from_VCF.py <input_vcf> <subtype> <sample_id> <transcripts_file> <token> <output_dir> <output>')
	sys.exit()
	
chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

if os.path.isfile(input_vcf):
	if input_vcf.endswith('.gz'):
		opened_input_vcf = gzip.open(input_vcf, 'rt', encoding='utf-8')
	else:
		opened_input_vcf = open(input_vcf, 'rt', encoding='utf-8')
else:
	print('cannot open '+input_vcf)
	sys.exit()
	
if os.path.isfile(transcripts_file):
	if transcripts_file.endswith('.gz'):
		opened_transcripts_file = gzip.open(transcripts_file, 'rt', encoding='utf-8')
	else:
		opened_transcripts_file = open(transcripts_file, 'rt', encoding='utf-8')
else:
	print('cannot open '+transcripts_file)
	sys.exit()
	
if os.path.isfile(token):
	if token.endswith('.gz'):
		opened_token = gzip.open(token, 'rt', encoding='utf-8')
	else:
		opened_token = open(token, 'rt', encoding='utf-8')
else:
	print('cannot open '+token)
	sys.exit()
	
token_str = str(opened_token.readline().strip())

print('found token '+token_str)

transcript_list = []

for line in opened_transcripts_file:
	line = line.strip()
	fields = line.split('\t')
	transcript = fields[0]
	
	if transcript not in transcript_list:
		transcript_list.append(transcript)
		
# list of nonsynonymous ensembl consequence terms - "promoter_variant" added to account for any TERT promoter mutations (which should be treated as nonsynonymous driver mutations)
nonsynonymous_list = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification',
'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification', 'inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant', 'promoter_variant']

# broaden list
consequence_broad_list = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification',
'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification', 'inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant', 'promoter_variant',
'splice_region_variant', 'splice_donor_5th_base_variant', 'splice_donor_region_variant', 'splice_polypyrimidine_tract_variant', 'incomplete_terminal_codon_variant', 'start_retained_variant', 'stop_retained_variant', 'synonymous_variant', 
'coding_sequence_variant', 'mature_miRNA_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', 'feature_truncation']

# rank consequences by severity according to https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
consequence_severity_dict = {"transcript_ablation":1, "splice_acceptor_variant":2, "splice_donor_variant":3, "stop_gained":4, "frameshift_variant":5, "stop_lost":6, "start_lost":7, "transcript_amplification":8,
"inframe_insertion":9, "inframe_deletion":10, "missense_variant":11, "protein_altering_variant":12, "splice_region_variant":13, "splice_donor_5th_base_variant":14, "splice_donor_region_variant":15,
"splice_pyrimidine_tract_variant":16, "incomplete_terminal_codon_variant":17, "start_retained_variant":18, "stop_retained_variant":19, "synonymous_variant":20, "coding_sequence_variant":21, 
"mature_miRNA_variant":22, "5_prime_UTR_variant":23, "3_prime_UTR_variant":24, "non_coding_transcript_exon_variant":25, "intron_variant":26, "NMD_transcript_variant":27, "non_coding_transcript_variant":28,
"upstream_gene_variant":29, "downstream_gene_variant":30, "TFBS_ablation":31, "TFBS_amplification":32, "TF_binding_site_variant":33, "regulatory_region_ablation":34, "regulatory_region_amplification":35,
"feature_elongation":36, "regulatory_region_variant":37, "feature_truncation":38, "intergenic_variant":39
}

oncokbapibearertoken = token_str

annotation = []

def makeoncokbgetrequest(url):
	headers = {
		'Content-Type': 'application/json',
		'Authorization': 'Bearer %s' % oncokbapibearertoken
	}
	return requests.get(url, headers=headers)
	
# write to output

opened_output = open(output, 'w')
output_writer = csv.writer(opened_output, delimiter = '\t')
output_writer.writerow(['gene', 'transcript', 'HGVSg', 'HGVSc', 'HGVSp', 'consequence'])

gene_mutations_list_dict = {}
mutation_count_dict = {}
transcript_dict = {}
protein_pos_dict = {}
label_dict = {}
carrier_dict = {}
gene_count_dict = {}
consequence_dict = {}

vep_lookup_dict = {}

for line in opened_input_vcf:
	# lookup VEP annotation orderings in VCF header
	if line.startswith('##INFO=<ID=CSQ,Number='):
		temp1 = line.split(':')
		temp2 = temp1[1].strip().split('"')
		vep_fields = temp2[0].strip().split('|')
		
		vep_ticker = 0
		
		for anno in vep_fields:
			vep_lookup_dict[anno] = vep_ticker
			vep_ticker += 1
			
	if line.startswith('##') or line.startswith('#CHROM'):
		continue
		
	line = line.strip()
	fields = line.split('\t')
	chr, pos, id, ref, alt, qual, filter, info, format = fields.split('\t', 8)
		
	# only keep PASS mutations in chr1-22, X, Y
	if chr not in chromosomes or (filter != "PASS" and filter != "NA"):
		continue
			
	# begin to parse VEP INFO fields for variant annotation
	info_split = info.split(';')
		
	for query in info_split:
		if query.startswith('CSQ='):
			vep_split = query.split(',')
				
	transcript_chosen = 'NA'
	consequence_chosen = 'NA'
	IMPACT_chosen = 'NA'
	Gene_chosen = 'NA'
	SYMBOL_chosen = 'NA'
	HGVSc_chosen = 'NA'
	HGVSp_chosen = 'NA'
	HGVSg_chosen = 'NA'
		
	# returns most damaging consequence in query transcript
	for query in vep_split:
		vep_string_split = query.split('|')
			
		canonical_flag = vep_string_split[vep_lookup_dict['CANONICAL']]
		consequence = vep_string_split[vep_lookup_dict['Consequence']]
		IMPACT = vep_string_split[vep_lookup_dict['IMPACT']]
		Gene = vep_string_split[vep_lookup_dict['Gene']]
		Transcript = vep_string_split[vep_lookup_dict['Feature']]
		SYMBOL = vep_string_split[vep_lookup_dict['SYMBOL']]
		HGVSc = vep_string_split[vep_lookup_dict['HGVSc']]
		HGVSp = vep_string_split[vep_lookup_dict['HGVSp']]
		HGVSg = vep_string_split[vep_lookup_dict['HGVSg']]
			
		consequence_split = consequence.split('&')
			
		if Transcript in transcript_list:
			for consequence in consequence_split:
				if consequence_chosen == 'NA':
					consequence_chosen = consequence
					transcript_chosen = Transcript
					IMPACT_chosen = IMPACT
					Gene_chosen = Gene
					SYMBOL_chosen = SYMBOL
					HGVSc_chosen = HGVSc
					HGVSp_chosen = HGVSp
					HGVSg_chosen = HGVSg
							
				elif int(consequence_severity_dict[consequence]) < int(consequence_severity_dict[consequence_chosen]):
					consequence_chosen = consequence
					transcript_chosen = transcript
					IMPACT_chosen = IMPACT
					Gene_chosen = Gene
					SYMBOL_chosen = SYMBOL
					HGVSc_chosen = HGVSc
					HGVSp_chosen = HGVSp
					HGVSg_chosen = HGVSg
							
				if HGVSp_chosen == '':
					HGVSp_chosen = 'NA'
							
		# rescuing TERT promoter mutations
		if chr == 'chr5' and (pos == '1295113' or pos == '1295135'):
			SYMBOL_chosen = 'TERT'
			transcript_chosen = 'ENST00000310581'
			consequence_chosen = 'promoter_variant'
			HGVSp_chosen = 'NA'
			HGVSc_chosen = 'NA'
					
		if ref in ['A', 'C', 'G', 'T'] and alt in ['A', 'C', 'G', 'T']:
			alt_type = 'SNV'
		else:
			alt_type = 'INDEL'
				
		# output nonsynonymous mutations in query genes
		if transcript_chosen in transcript_list and consequence_chosen in nonsynonymous_list:
				
			output_writer.writerow([str(SYMBOL_chosen), str(transcript_chosen), str(chr), str(pos), str(ref), str(alt), str(HGVSg_chosen), str(HGVSc_chosen), str(HGVSp_chosen), str(consequence_chose)])
				
			# query OncoKB API to get annotation
				
			json_output = output_dir+'/'+str(chr)+'_'+str(pos)+'_'+str(ref)+'_'+str(alt)+'_'+str(SYMBOL_chosen)+'_'+str(HGVSp_chosen)+'_oncoKB_query.json'
				
			if consequence_chosen in nonsynonymous_list and HGVSp_chosen != "NA":
				HGVSp_split = HGVSp_chosen.split(':')
				HGVSp_shortened = HGVSp_split[1]
					
				protein_geturl = 'https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?'
				protein_geturl += 'hugoSymbol='+str(SYMBOL_chosen)+'&alteration='+str(HGVSp_shortened)+'&consequence='+str(consequence_chosen)
				protein_geturl += '&referenceGenome=GRCh38&tumorType='+str(subtype)
					
			elif consequence_chosen in nonsynonymous_list:
				protein_geturl = 'https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?'
				protein_geturl += 'hugoSymbol='+str(SYMBOL_chosen)+'&alteration='+str(HGVSp_chosen)+'&consequence='+str(consequence_chosen)
				protein_geturl += '&referenceGenome=GRCh38&tumorType='+str(subtype)
					
			if alt_type == 'INDEL':
				geturl = 'https://www.oncokb.org/api/v1/annotate/mutations/byHGVSg?'
				geturl += 'hgvsg='+str(HGVSg_chosen)
				geturl += '&referenceGenome=GRCh38&tumorType='+str(subtype)
					
			else:
				geturl = 'https://www.oncokb.org/api/v1/annotate/mutations/byHGVSg?'
				geturl += 'hgvsg='+str(chr)+':g.'+str(pos)+str(ref)+'%3E'+str(alt)
				geturl += '&referenceGenome=GRCh38&tumorType='+str(subtype)
				
			query_ticker = 1
			query_success = 'NO'
			knownEffect = ''
			hugoSymbol = ''
				
			for i in range (1, 11):
				if i == 10 and query_success != "YES":
					print('failed 10 times so exiting, perhaps try again later')
					sys.exit()
						
				elif query_success == 'YES':
					continue
						
				elif i > 1 and getresponse.status_code == 200 and knownEffect != '' and hugoSymbol != '':
					annotation.append(getresponse.json())
					query_success = 'YES'
					continue
						
				else:
					if len(geturl) > 200:
						getresponse = makeoncokbgetrequest(protein_geturl)
					else:
						getresponse = makeoncokbgetrequest(geturl)
							
					time.sleep(1)
						
					if getresponse.status_code != 200:
						getresponse = makeoncokbgetrequest(protein_geturl)
							
					mutation_dict = getresponse.json()
					knownEffect = mutation_dict['mutationEffect']['knownEffect']
					hugoSymbol = mutation_dict['query']['hugoSymbol']
						
combined_json = output_dir+'/'+plate_tumour+'_combined_mutations_OncoKB_query.json'
opened_combined_json = open(combined_json, 'w')

for anno in annotation:
	json.dump(anno, opened_combined_json)
	opened_combined_json.write('\n')
	
opened_combined_json.close()
opened_output.close()

