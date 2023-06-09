#!/usr/bin/env python3

### generate file containing interpro domain information suitable for plotting (for use with MAFtools etc)
### Ben Kinnersley (b.kinnersley@ucl.ac.uk)
### Command line: python make_interpro_file_for_plotting.py <regions_interpro> <interpro_entry> <protein_fasta> <output>

import os
import sys
import gzip
import csv

if len(sys.argv) == 5:
	regions_interpro = sys.argv[1]
	interpro_entry = sys.argv[2]
	protein_fasta = sys.argv[3]
	output = sys.argv[4]
else:
	print('insufficient arguments provided - python make_interpro_file_for_plotting.py <regions_interpro> <interpro_entry> <protein_fasta> <output>')
	sys.exit()
	
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']

if os.path.isfile(regions_interpro):
	if regions_interpro.endswith('.gz'):
  		opened_regions_interpro = gzip.open(regions_interpro,'rt', encoding='utf-8')
	else:
		opened_regions_interpro = open(regions_interpro, 'rt', encoding='utf-8')
	print('reading interpro regions '+regions_interpro)
else:
	print('cannot open '+regions_interpro) 
	sys.exit()
	
if os.path.isfile(interpro_entry):
	if interpro_entry.endswith('.gz'):
  		opened_interpro_entry = gzip.open(interpro_entry,'rt', encoding='utf-8')
	else:
		opened_interpro_entry = open(interpro_entry, 'rt', encoding='utf-8')
	print('reading interpro entry file '+interpro_entry)
else:
	print('cannot open '+interpro_entry) 
	sys.exit()

if os.path.isfile(protein_fasta):
	if protein_fasta.endswith('.gz'):
  		opened_protein_fasta = gzip.open(protein_fasta,'rt', encoding='utf-8')
	else:
		opened_protein_fasta = open(protein_fasta, 'rt', encoding='utf-8')
	print('reading protein fasta '+protein_fasta)
else:
	print('cannot open '+protein_fasta) 
	sys.exit()
		
opened_output = open(output, 'w')
output_writer = csv.writer(opened_output, delimiter = '\t')

print('writing to output '+output)

canonical_transcripts_list = []
all_transcripts_list = []
gene_symbol_dict = {}

opened_regions_interpro.readline()

interpro_domain_list = []
interpro_per_transcript_dict = {}
interpro_desc_short_dict = {}
interpro_desc_full_dict = {}

for line in opened_regions_interpro:
	gene_id, transcript_id, canonical, protein_id, interpro_id, interpro_desc_short, interpro_desc_full, interpro_start, interpro_end, gene_symbol = line.split('\t', 9)
	gene_symbol = gene_symbol.strip()
	
	# get list of canonical transcripts
	if canonical == '1':
		if transcript_id not in canonical_transcripts_list:
			canonical_transcripts_list.append(transcript_id)
	
	if transcript_id not in all_transcripts_list:
		all_transcripts_list.append(transcript_id)
	
	if interpro_id not in interpro_domain_list and interpro_id != '':
		interpro_domain_list.append(interpro_id)
	
	if interpro_id != '':
		interpro_desc_short_dict[interpro_id] = interpro_desc_short
		interpro_desc_full_dict[interpro_id] = interpro_desc_full
	
	gene_symbol_dict[transcript_id] = gene_symbol
	
	if transcript_id not in interpro_per_transcript_dict and interpro_id != '':
		interpro_per_transcript_dict[transcript_id] = str(interpro_id)+':'+str(interpro_start)+':'+str(interpro_end)
	elif interpro_id != '':
		interpro_per_transcript_dict[transcript_id] = interpro_per_transcript_dict[transcript_id]+'/'+str(interpro_id)+':'+str(interpro_start)+':'+str(interpro_end)

interpro_type_dict = {}		
		
for line in opened_interpro_entry:
	line = line.strip()
	interpro_id, interpro_type, interpro_desc = line.split('\t', 2)
	
	interpro_type_dict[interpro_id] = interpro_type
		
protein_sequence_dict = {}
protein_name_dict = {}

for line in opened_protein_fasta:
	if line.startswith('>'):
		line_split = line.split(' ')
		
		for query in line_split:
			if query.startswith('>'):
				query_split = query.split('>')
				protein_split = query_split[1].split('.')
				protein = protein_split[0]
				
			if query.startswith('chromosome:GRCh38') and not query.startswith('chromosomes'):
				query_split = query.split(':')
				chromosome = query_split[2]
				
			if query.startswith('transcript:'):
				query_split = query.split(':')
				transcript_split = query_split[1].split('.')
				transcript = transcript_split[0]
				
		amino_acid_string = ''
	else:
		amino_acid_string = amino_acid_string+line.rstrip()
		
		if transcript in all_transcripts_list and chromosome in chromosomes:
			protein_name_dict[transcript] = protein
			protein_sequence_dict[transcript] = amino_acid_string
			
output_writer.writerow(['HGNC', 'refseq.ID', 'protein.ID', 'aa.length', 'Start', 'End', 'pfam', 'Label', 'Description', 'Type'])

for transcript in all_transcripts_list:
	if transcript in interpro_per_transcript_dict:
		query_split = interpro_per_transcript_dict[transcript].split('/')
		for query in query_split:
			interpro_split = query.split(':')
			interpro_domain = interpro_split[0]
			interpro_start = interpro_split[1]
			interpro_end = interpro_split[2]
				
			output_writer.writerow([str(gene_symbol_dict[transcript]), str(transcript), str(protein_name_dict[transcript]), str(len(protein_sequence_dict[transcript])),
				str(interpro_start), str(interpro_end), str(interpro_domain), str(interpro_desc_short_dict[interpro_domain]), str(interpro_desc_full_dict[interpro_domain]),
				str(interpro_type_dict[interpro_domain])])
			
	elif transcript in protein_sequence_dict:
		output_writer.writerow([str(gene_symbol_dict[transcript]), str(transcript), str(protein_name_dict[transcript]), str(len(protein_sequence_dict[transcript])),
			'NA', 'NA', 'NA', 'NA', 'NA', 'NA'])
