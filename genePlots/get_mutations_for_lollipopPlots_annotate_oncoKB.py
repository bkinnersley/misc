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
	transcripts = fields[0].strip()
	
	if transcript not in oncokb_transcript_list:
		oncokb_transcript_list.append(transcript)
		
	if transcript not in transcript_list:
		print(str(transcript)+' (oncoKB transcript) is not canonical'

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
				transcript = transcript_split[0]
		
	
	
	
