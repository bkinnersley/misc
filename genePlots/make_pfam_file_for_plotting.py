#!/usr/bin/env python3

### generate file containing pfam domain information suitable for plotting (for use with MAFtools etc)
### Ben Kinnersley (b.kinnersley@ucl.ac.uk)
### Command line: python make_pfam_file_for_plotting.py <regions_pfam> <pfam_names> <protein_fasta> <pfam_dead> <output>

import os
import sys
import gzip

if len(sys.argv) == 6:
	regions_pfam = sys.argv[1]
	pfam_names = sys.argv[2]
	protein_fasta = sys.argv[3]
	pfam_dead = sys.argv[4]
	output = sys.argv[5]
else:
	print('insufficient arguments provided - python make_pfam_file_for_plotting.py <regions_pfam> <pfam_names> <protein_fasta> <pfam_dead> <output>')
	sys.exit()

chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']

if os.path.isfile(regions_pfam):
	if regions_pfam.endswith('.gz'):
  		opened_regions_pfam = gzip.open(regions_pfam,'rb')
	else:
		opened_regions_pfam = gzip.open(regions_pfam, 'rb')
	print('reading pfam regions '+regions_pfam)
else:
	print('cannot open '+regions_pfam) 
	sys.exit()
	
if os.path.isfile(pfam_names):
	if pfam_names.endswith('.gz'):
  		opened_pfam_names = gzip.open(pfam_names,'rb')
	else:
		opened_pfam_names = gzip.open(pfam_names, 'rb')
	print('reading pfam names '+pfam_names)
else:
	print('cannot open '+pfam_names) 
	sys.exit()

if os.path.isfile(protein_fasta):
	if protein_fasta.endswith('.gz'):
  		opened_protein_fasta = gzip.open(protein_fasta,'rb')
	else:
		opened_protein_fasta = gzip.open(protein_fasta, 'rb')
	print('reading protein fasta '+protein_fasta)
else:
	print('cannot open '+protein_fasta) 
	sys.exit()
	
if os.path.isfile(pfam_dead):
	if pfam_dead.endswith('.gz'):
  		opened_pfam_dead = gzip.open(pfam_dead,'rb')
	else:
		opened_pfam_names = gzip.open(pfam_dead, 'rb')
	print('reading dead pfam information '+pfam_dead)
else:
	print('cannot open '+pfam_dead) 
	sys.exit()
	
opened_output = open(output, 'w')
output_writer = csv.writer(opened_output, delimiter = '\t')
#output_writer.writerow(['gene', 'transcript', 'HGVSg', 'HGVSc', 'HGVSp', 'consequence'])

print('writing to output '+output)

canonical_transcripts_list = []
all_transcripts_list = []
gene_symbol_dict = {}

opened_regions_pfam.readline()

pfam_domain_list = []
pfam_per_transcript_dict = {}

for line in opened_regions_pfam:
	gene_id, transcript_id, pfam_id, pfam_start, pfam_end, protein_id, canonical, gene_symbol = line.split('\t', 7)
	
	# get list of canonical transcripts
	if int(canonical) == 1:
		if transcript_id not in canonical_transcripts_list:
			canonical_transcripts_list.append(transcript_id)
	
	if transcript_id not in all_transcripts_list:
		all_transcripts_list.append(transcript_id)
	
	gene_symbol_dict[transcript_id] = gene_symbol
	
	if transcript_id not in pfam_per_transcript_dict:
		pfam_per_transcript_dict[transcript_id] = str(pfam_id)+':'+str(pfam_start)+':'+str(pfam_end)
	else:
		pfam_per_transcript_dict[transcript_id] = pfam_per_transcript_dict[transcript_id]+'/'+str(pfam_id)+':'+str(pfam_start)+':'+str(pfam_end)
	
# lookup deprecated/changed pfam domain annotations

pfam_start_id = ''
pfam_replacement_id = ''
pfam_replacement_dict = {}
pfam_dead_list = []

for line in opened_pfam_dead:
	if line.startswith('#=GF AC'):
		fields = line.split()
		pfam_start_id = fields[2].strip()
		
	if line.startswith('#=GF FW'):
		fields = line.split()
		
		if len(fields) == 3:
			pfam_replacement_id = fields[2].strip()
			pfam_replacement_dict[pfam_replacement_id] = pfam_start_id
			
		else:
			pfam_dead_list.append(pfam_start_id)
			pfam_replacement_id = 'MISSING'
	if line.startswith('//'):
		pfam_start_id = ''
		pfam_replacement_id = ''
		
pfam_family_dict = {}
pfam_description_dict = {}

for line in opened_pfam_names:
	fields = line.split('\t')
	pfam_domain = fields[0].strip()
	pfam_family = fields[3].strip()
	pfam_description = fields[4].strip()
	
	if pfam_domain in pfam_replacement_dict:
		pfam_family_dict[pfam_replacement_dict[pfam_domain]] = pfam_family
		pfam_description_dict[pfam_replacement_dict[pfam_domain]] = pfam_description
		
	if pfam_domain in pfam_domain_list:
		pfam_family_dict[pfam_domain] = pfam_family
		pfam_description_dict[pfam_domain] = pfam_description
		
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
		
		if transcript in canonical_transcripts_list and chromosome in chromosomes:
			protein_name_dict[transcript] = protein
			protein_sequence_dict[transcript] = amino_acid_string
			
output_writer.writerow(['HGNC', 'refseq.ID', 'protein.ID', 'aa.length', 'Start', 'End', 'pfam', 'Label', 'Description'])

for transcript in all_transcripts_list:
	if transcript in pfam_per_transcript_dict:
		query_split = pfam_per_transcript_dict[transcript].split('/')
		for query in query_split:
			pfam_split = query.split(':')
			pfam_domain = pfam_split[0]
			pfam_start = pfam_split[1]
			pfam_end = pfam_split[2]
			
			if pfam_domain not in pfam_family_dict:
				pfam_family_dict[pfam_domain] = str(pfam_domain)+'_NA'
				pfam_description_dict[pfam_domain] = str(pfam_domain)+'_NA'
				
			output_writer.writerow([gene_symbol_dict[transcript], str(transcript), str(protein_name_dict[transcript]), str(len(protein_sequence_dict[transcript])),
				str(pfam_start), str(pfam_end), str(pfam_domain), str(pfam_family_dict[pfam_domain]), str(pfam_description_dict[pfam_domain])])
			
		elif transcript in protein_sequence_dict:
			output_write.writerow([gene_symbol_dict[transcript], str(transcript), str(protein_name_dict[transcript]), str(len(protein_sequence_dict[transcript])),
				'NA', 'NA', 'NA', 'NA', 'NA'])
									
	
