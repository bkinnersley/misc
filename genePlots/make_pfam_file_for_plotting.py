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
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
