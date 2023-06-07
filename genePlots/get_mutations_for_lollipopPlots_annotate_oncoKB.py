#!/usr/bin/env python3

### extract nonsynonymous mutations from user supplied transcripts (e.g. ensembl canonical transcripts) from a query set of genes and tumour samples and add oncokb annotation
### Ben Kinnersley (b.kinnersley@ucl.ac.uk)
### Command line: python get_mutations_for_lollipopPlots_annotate_oncoKB.py <input_samples> <vcf_dir> <gene_list> <pfam_file> <transcripts_file> <protein_fasta> <oncokb_folder> <oncokb_genes> <output>

import os
import sys
import gzip 
import re
import csv

if len(sys.argv) == 10:
	input_samples = sys.argv[1]
	vcf_dir = sys.argv[2]
	gene_list = sys.argv[3]
	pfam_file = sys.argv[4]
	transcripts_file = sys.argv[5]
	protein_fasta = sys.argv[6]
	oncokb_folder = sys.argv[7]
	oncokb_genes = sys.argv[8]
	output = sys.argv[9]
else:
	print('insufficient arguments provided - python get_mutations_for_lollipopPlots_annotate_oncoKB.py <input_samples> <vcf_dir> <gene_list> <pfam_file> <transcripts_file> <protein_fasta> <oncokb_folder> <oncokb_genes> <output>')
	sys.exit()

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

if os.path.isfile(input_samples):
	if input_samples.endswith('.gz'):
  		opened_input_samples = gzip.open(input_samples,'rt', encoding='utf-8')
	else:
		opened_input_samples = open(input_samples, 'rt', encoding='utf-8')
	print('reading samples '+input_samples)
else:
	print('cannot open '+input_samples) 
	sys.exit()
  
if os.path.isfile(gene_list):
	if gene_list.endswith('.gz'):
  		opened_gene_list = gzip.open(gene_list,'rt', encoding='utf-8')
	else:
		opened_gene_list = open(gene_list, 'rt', encoding='utf-8')
	print('reading genes '+gene_list)
else:
	print('cannot open '+gene_list) 
	sys.exit()
	
if os.path.isfile(pfam_file):
	if pfam_file.endswith('.gz'):
  		opened_pfam_file = gzip.open(pfam_file,'rt', encoding='utf-8')
	else:
		opened_pfam_file = open(pfam_file, 'rt', encoding='utf-8')
	print('reading pfam information '+pfam_file)
else:
	print('cannot open '+pfam_file) 
	sys.exit()

if os.path.isfile(transcripts_file):
	if transcripts_file.endswith('.gz'):
  		opened_transcripts_file = gzip.open(transcripts_file,'rt', encoding='utf-8')
	else:
		opened_transcripts_file = open(transcripts_file, 'rt', encoding='utf-8')
	print('reading transcripts '+transcripts_file)
else:
	print('cannot open '+transcripts_file) 
	sys.exit()
  
if os.path.isfile(protein_fasta):
	if protein_fasta.endswith('.gz'):
  		opened_protein_fasta = gzip.open(protein_fasta,'rt', encoding='utf-8')
	else:
		opened_protein_fasta = open(protein_fasta, 'rt', encoding='utf-8')
	print('reading fasta information '+protein_fasta)
else:
	print('cannot open '+protein_fasta)
	
if os.path.isfile(oncokb_genes):
	if oncokb_genes.endswith('.gz'):
  		opened_oncokb_genes = gzip.open(oncokb_genes,'rt', encoding='utf-8')
	else:
		opened_oncokb_genes = open(oncokb_genes, 'rt', encoding='utf-8')
	print('reading oncokb gene information '+oncokb_genes)
else:
	print('cannot open '+oncokb_genes)	

	sys.exit()
  
if os.path.isfile(pfam_file):
	if pfam_file.endswith('.gz'):
  		opened_pfam_file = gzip.open(pfam_file,'rt', encoding='utf-8')
	else:
		opened_pfam_file = open(pfam_file, 'rt', encoding='utf-8')
	print('reading pfam information '+pfam_file)
else:
	print('cannot open '+pfam_file) 
	sys.exit()
	
query_gene_list = []

for line in opened_gene_list:
	fields = line.split('\t')
	gene = fields[0].strip()
	
	if gene not in query_gene_list:
		query_gene_list.append(gene)

sample_count = 0

transcript_dict = {}
ensembl_gene_list = []
transcript_list = []

for line in opened_transcripts_file:
	ensembl_gene_id, ensembl_transcript_id, gene_symbol = line.split('\t', 2)
	gene_symbol = gene_symbol.strip()
	
	if gene_symbol in query_gene_list:
		ensembl_gene_list.append(ensembl_gene_id)
	
	if ensembl_transcript_id not in transcript_list:
		transcript_list.append(ensembl_transcript_id)
		
	transcript_dict[ensembl_gene_id] = ensembl_transcript_id

oncokb_transcript_list = []
	
for line in opened_oncokb_genes:
	fields = line.split('\t')
	transcript = fields[0].strip()
	
	if transcript not in oncokb_transcript_list:
		oncokb_transcript_list.append(transcript)
		
	if transcript not in transcript_list:
		print(str(transcript)+' (oncoKB transcript) is not canonical')

protein_sequence_dict = {}
protein_name_dict = {}
protein_size_dict = {}
		
for line in opened_protein_fasta:
	if line.startswith('>'):
		line_split = line.split(' ')
		
		for query in line_split:
			if query.startswith('>'):
				query_split = query.split('>')
				protein_split = query_split[1].split('.')
				protein = protein_split[0]
			
			if query.startswith('transcript:'):
				query_split = query.split(':')
				transcript_split = query_split[1].split('.')
				transcript = transcript_split[0]  
		    
		amino_acid_string = ''
	else:
		amino_acid_string = amino_acid_string+line.rstrip()
		
		if transcript in oncokb_transcript_list:
			  protein_name_dict[transcript] = protein
			  protein_sequence_dict[transcript] = amino_acid_string
			  protein_size_dict[transcript] = len(amino_acid_string)
			  
opened_pfam_file.readline()

for line in opened_pfam_file:
	line = line.strip()
	gene_symbol, transcript, protein, protein_length = line.split('\t', 3)
	
	protein_size_dict[transcript] = protein_length
			  
# list of nonsynonymous ensembl consequence terms - "promoter_variant" added to account for any TERT promoter mutations (which should be treated as nonsynonymous driver mutations)
nonsynonymous_list = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification',
'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification', 'inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant', 'promoter_variant']

# rank consequences by severity according to https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
consequence_severity_dict = {"transcript_ablation":1, "splice_acceptor_variant":2, "splice_donor_variant":3, "stop_gained":4, "frameshift_variant":5, "stop_lost":6, "start_lost":7, "transcript_amplification":8,
"inframe_insertion":9, "inframe_deletion":10, "missense_variant":11, "protein_altering_variant":12, "splice_region_variant":13, "splice_donor_5th_base_variant":14, "splice_donor_region_variant":15,
"splice_polypyrimidine_tract_variant":16, "incomplete_terminal_codon_variant":17, "start_retained_variant":18, "stop_retained_variant":19, "synonymous_variant":20, "coding_sequence_variant":21, 
"mature_miRNA_variant":22, "5_prime_UTR_variant":23, "3_prime_UTR_variant":24, "non_coding_transcript_exon_variant":25, "intron_variant":26, "NMD_transcript_variant":27, "non_coding_transcript_variant":28,
"upstream_gene_variant":29, "downstream_gene_variant":30, "TFBS_ablation":31, "TFBS_amplification":32, "TF_binding_site_variant":33, "regulatory_region_ablation":34, "regulatory_region_amplification":35,
"feature_elongation":36, "regulatory_region_variant":37, "feature_truncation":38, "intergenic_variant":39
}

gene_mutations_list_dict = {}
mutation_count_dict = {}
transcript_dict = {}
protein_pos_dict = {}
label_dict = {}
carrier_dict = {}
gene_count_dict = {}
consequence_dict = {}
			  
oncogenic_dict = {}
knownEffect_dict = {}
hotspot_dict = {}
vus_dict = {}
oncogenic_output_dict = {}
knownEffect_output_dict = {}
hotspot_output_dict = {}
vus_output_dict = {}

for line in opened_input_samples:
	line = line.strip()
	fields = line.split('\t')
	sample_id = fields[0]
	
	# to do: improve parsing of paths
	input_oncokb = oncokb_folder+'/'+str(sample_id)+'_combined_mutations_OncoKB_query_parsed.tsv'
	
	if os.path.isfile(input_oncokb):
		if input_oncokb.endswith('.gz'):
			  opened_input_oncokb = gzip.open(input_oncokb, 'rt')
		else:
			  opened_input_oncokb = open(input_oncokb, 'rt')
	else:
		print('could not open file ')+input_oncokb
		continue
	
	opened_input_oncokb.readline()
	
	for line in opened_input_oncokb:
		line = line.strip()
		
		fields = line.split('\t')
		
		oncogenic = fields[5]
		knownEffect = fields[6]
		hugoSymbol = fields[12]
		alteration = fields[13]
		hgvs = fields[15]
		consequence = fields[16]
		hotspot = fields[18]
		vus = fields[31]
		
		# for SNVs match by HGVSg, for nonsynonymous indels match by gene and protein change (to do: ensure indels have proper HGVSg matching)
		if hgvs != "None":
			oncogenic_dict[hgvs] = oncogenic
			knownEffect_dict[hgvs] = knownEffect
			hotspot_dict[hgvs] = hotspot
			vus_dict[hgvs] = vus
		else:
			query = hugoSymbol+'_p.'+alteration+'_'+consequence
			
			oncogenic_dict[query] = oncogenic
			knownEffect_dict[query] = oncogenic
			hotspot_dict[query] = hotspot
			vus_dict[query] = vus
			
	input_vcf = vcf_dir+'/GCGR-'+str(sample_id)+'/'+str(sample_id)+'_filt_pass.VEP.vcf.gz'
	
	if os.path.isfile(input_vcf):
		if input_vcf.endswith('.gz'):
			opened_input_vcf = gzip.open(input_vcf, 'rt')
		else:
			opened_input_vcf = open(input_vcf, 'rt')
		
		sample_count += 1
		print(str(sample_count)+': reading vcf '+input_vcf)
	else:
		print('could not read file '+vcf)
		continue
	
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
		chr, pos, id, ref, alt, qual, filter, info, format = line.split('\t', 8)
			
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
		if Gene_chosen in ensembl_gene_list and consequence_chosen in nonsynonymous_list:
			if Gene_chosen in transcript_dict:
				transcript = transcript_dict[Gene_chosen]
				
			mutation = chr+'_'+pos+'_'+ref+'_'+alt
			
			# if mutation is protein-coding take protein position from HGVSp annotation, otherwise approximate for splice donor/acceptor mutations using HGVSc change
			
			if HGVSp_chosen != '' and HGVSp_chosen != 'NA':
				HGVSp_split = HGVSp_chosen.split(':')
				HGVSp_shortened = HGVSp_split[1]
				
				pattern = r"[p\.]+[a-zA-Z]+([0-9]+)[\_]*[a-zA-Z\?\=\%]+"
				
				result = re.search(pattern, HGVSp_shortened)
				
				if result.group() is not None:
					pass
				else:
					print('check parsing of '+str(HGVSp_chosen)+' in mutation '+str(mutation))
					sys.exit()
				
				protein_pos = result.group(1)
				
				HGVSc_split = HGVSc_chosen.split('.')
				transcript = HGVSc_split[0]
				
				label = HGVSp_shortened
				
			elif HGVSc_chosen != '' and HGVSc_chosen != 'NA':
				HGVSc_split = HGVSc_chosen.split('.')
				transcript = HGVSc_split[0]
				
				# if 5' UTR splice mutation set to protein position 1
				if HGVSc_split[2].startswith('-'):
					protein_pos = 1
				# if 3' UTR splice mutation set to last protein position
				elif HGVSc_split[2].startswith('*'):
					protein_pos = protein_size_dict[transcript]
				elif '-' in HGVSc_split[2]:
					HGVSc_temp = HGVSc_split[2].split('-')
					pos_split = HGVSc_temp[0].split('_')
					transcript_pos = pos_split[0]
					
					if '+' in pos_split[0]:
						pos_split_again = pos_split[0].split('+')
						transcript_pos = pos_split_again[0]
					else:
						transcript_pos = pos_split[0]
						
					protein_pos = int(round(float(transcript_pos)/3))
				elif '+' in HGVSc_split[2]:
					HGVSc_temp = HGVSc_split[2].split('+')
					pos_split = HGVSc_temp[0].split('_')
					transcript_pos = pos_split[0]
					
					protein_pos = int(round(float(transcript_pos)/3))
				else:
					'check! '+HGVSc_chosen
				
				label = 'c.'+HGVSc_split[2]
			
			elif SYMBOL_chosen == 'TERT' and consequence_chosen == 'promoter_variant':
				print('seen TERT promoter variant '+str(mutation))
				protein_pos = 1
				
				if pos == '1295113':
					label = 'C250T'
				elif pos == '1295135':
					label = 'C228T'
			else:
				print(input_vcf)
				print('ignoring query mutation with missing HGVSc and HGVSp '+str(mutation)+' ('+str(SYMBOL_chosen)+') - '+str(consequence_chosen))
				continue
			
			protein_pos_dict[mutation] = protein_pos
			label_dict[mutation] = label
			transcript_dict[mutation] = transcript
			consequence_dict[mutation] = consequence_chosen
			
			if SYMBOL_chosen not in gene_mutations_list_dict:
				gene_mutations_list_dict[SYMBOL_chosen] = []
			
			if mutation not in gene_mutations_list_dict[SYMBOL_chosen]:
				gene_mutations_list_dict[SYMBOL_chosen].append(mutation)
				
			if mutation not in carrier_dict:
				carrier_dict[mutation] = sample_id
			else:
				carrier_dict[mutation] = carrier_dict[mutation]+';'+sample_id
				
			if mutation in mutation_count_dict:
				mutation_count_dict[mutation] += 1
			else:
				mutation_count_dict[mutation] = 1
				
			if SYMBOL_chosen in gene_count_dict:
				gene_count_dict[SYMBOL_chosen] += 1
			else:
				gene_count_dict[SYMBOL_chosen] = 1
			  
			# try to link back with oncoKB mutations
			HGVSp_mut_query = SYMBOL_chosen+'_'+label+'_'+consequence_chosen
			
			if HGVSg_chosen in oncogenic_dict:
				oncogenic_output_dict[mutation] = oncogenic_dict[HGVSg_chosen]
				knownEffect_output_dict[mutation] = knownEffect_dict[HGVSg_chosen]
				hotspot_output_dict[mutation] = hotspot_dict[HGVSg_chosen]
				vus_output_dict[mutation] = vus_dict[HGVSg_chosen]
			
			elif HGVSp_mut_query in oncogenic_dict:
				oncogenic_output_dict[mutation] = oncogenic_dict[HGVSp_mut_query]
				knownEffect_output_dict[mutation] = knownEffect_dict[HGVSp_mut_query]
				hotspot_output_dict[mutation] = hotspot_dict[HGVSp_mut_query]
				vus_output_dict[mutation] = vus_dict[HGVSp_mut_query]
			else:
				oncogenic_output_dict[mutation] = 'NA'
				knownEffect_output_dict[mutation] = 'NA'
				hotspot_output_dict[mutation] = 'NA'
				vus_output_dict[mutation] = 'NA'

# writing to output

opened_output = open(output, 'w')
output_writer = csv.writer(opened_output, delimiter = '\t')
output_writer.writerow(['gene', 'transcript', 'mutation', 'count', 'protein_pos', 'label', 'consequence', 'oncogenic', 'knownEffect', 'hotspot', 'vus', 'carriers'])

for gene in query_gene_list:
	if gene not in gene_count_dict:
		print('0 mutations observed for '+gene)
		continue
	else:
		print(str(gene_count_dict[gene])+' mutations observed for '+gene)
		
	for mutation in gene_mutations_list_dict[gene]:
		output_writer.writerow([str(gene), str(transcript_dict[mutation]), str(mutation), str(mutation_count_dict[mutation]), str(protein_pos_dict[mutation]),
			str(label_dict[mutation]), str(consequence_dict[mutation]), str(oncogenic_output_dict[mutation]), str(knownEffect_output_dict[mutation]),
			str(hotspot_output_dict[mutation]), str(vus_output_dict[mutation]), str(carrier_dict[mutation])])
			
			
			  
			  
			  
