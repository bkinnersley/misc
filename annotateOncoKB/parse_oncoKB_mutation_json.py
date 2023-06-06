#!/usr/bin/env python3

### parse OncoKB mutation json ###
### Ben Kinnersley (b.kinnersley@ucl.ac.uk) ###
### Command line: python parse_oncoKB_mutation_json.py <input> <output>

import os
import sys
import gzip
import re
import requests
import json
import io
import csv

if len(sys.argv) == 3:
	input = sys.argv[1]
	output = sys.argv[2]
else:
	print('insufficient arguments provided - python parse_oncoKB_mutation_json.py <input> <output>')
	
if os.path.isfile(input):
	if input.endswith('.gz'):
		opened_input = gzip.open(input, encoding='utf-8')
	else:
		opened_input = io.open(input, encoding='utf-8')
	print('reading input: '+input)
else:
	print('cannot open '+input)
	sys.exit()
	
opened_output = open(output, 'w')
output_writer = csv.writer(opened_output, delimiter = '\t')
output_writer.writerow(['geneSummary', 'highestSensitiveLevel', 'dataVersion', 'tumorTypeSummary', 'highestPrognosticImplicationLevel', 'oncogenic', 'alleleExist', 'citations', 'knownEffect', 'proteinStart', 'proteinEnd', 'tumorType', 'hugoSymbol',
'alteration', 'entrezGeneId', 'hgvs', 'consequence', 'referenceGenome', 'hotspot', 'diagnosticImplications', 'prognosticSummary', 'lastUpdate', 'prognosticImplications', 'treatments', 'geneExist', 'variantExist', 'highestDiagnosticImplicationLevel',
'diagnosticSummary', 'otherSignificantSensitiveLevels', 'otherSignificantResistanceLevels', 'variantSummary', 'vus', 'highestResistanceLevel'])

for line in opened_input:
	line = line.decode('utf-8').strip()
	
	mutation_dict = json.loads(line)
	
	geneSummary = json.dumps(mutation_dict['geneSummary'], ensure_ascii = False)
	highestSensitiveLevel = mutation_dict['highestSensitiveLevel']
	dataVersion = mutation_dict['dataVersion']
	tumorTypeSummary = mutation_dict['tumorTypeSummary']
	highestPrognosticImplicationLevel = mutation_dict['highestPrognosticImplicationLevel']
	oncogenic = mutation_dict['oncogenic']
	alleleExist = mutation_dict['alleleExist']
	
	mutationEffect = mutation_dict['mutationEffect']
	citations = ','.join(map(str, mutation_dict['mutationEffect']['citations']['pmids']))
	knownEffect = mutation_dict['mutationEffect']['knownEffect']
	description = mutation_dict['mutationEffect']['description'] # seems to be empty
	
	query = mutation_dict['query']
	proteinStart = mutation_dict['query']['proteinStart']
	proteinEnd = mutation_dict['query']['proteinEnd']
	tumorType = mutation_dict['query']['tumorType']
	hugoSymbol = mutation_dict['query']['hugoSymbol']
	alteration = mutation_dict['query']['alteration']
	entrezGeneId = mutation_dict['query']['entrezGeneId']
	hgvs = mutation_dict['query']['hgvs']
	consequence = mutation_dict['query']['consequence']
	referenceGenome = mutation_dict['query']['referenceGenome']
	
	hotspot = mutation_dict['hotspot']
	diagnosticImplications = mutation_dict['diagnosticImplications']
	prognosticSummary = mutation_dict['prognosticSummary']
	lastUpdate = mutation_dict['lastUpdate']
	prognosticImplications = mutation_dict['prognosticImplications']
	treatments = json.dumps(mutation_dict['treatments'], ensure_ascii = False) # currently just dumping everything - TO DO: better parsing of treatment information
	geneExist = mutation_dict['geneExist']
	variantExist = mutation_dict['variantExist']
	highestDiagnosticImplicationLevel = mutation_dict['highestDiagnosticImplicationLevel']
	diagnosticSummary = mutation_dict['diagnosticSummary']
	otherSignificantSensitiveLevels = mutation_dict['otherSignificantSensitiveLevels']
	otherSignificantResistanceLevels = mutation_dict['otherSignificantResistanceLevels']
	variantSummary = mutation_dict['variantSummary']
	vus = mutation_dict['vus']
	highestResistanceLevel = mutation_dict['highestResistanceLevel']
	
	var_list = [geneSummary, highestSensitiveLevel, dataVersion, tumorTypeSummary, highestPrognosticImplicationLevel, oncogenic, alleleExist, citations, knownEffect, proteinStart, proteinEnd, tumorType, hugoSymbol, alteration,
	entrezGeneId, hgvs, consequence, referenceGenome, hotspot, diagnosticImplications, prognosticSummary, lastUpdate, prognosticImplications, treatments, geneExist, variantExist, highestDiagnosticImplicationLevel,
	diagnosticSummary, otherSignificantSensitiveLevels, otherSignificantResistanceLevels, variantSummary, vus, highestResistanceLevel]
	
	for var in var_list:
		if type(var) == unicode:
			var = var.encode('utf-8')
	
	#if type(geneSummary) == unicode:
	#	geneSummary = geneSummary.encode('utf-8')
		
	#if type(highestSensitiveLevel) == unicode:
	#	highestSensitiveLevel = highestSensitiveLevel.encode('utf-8')
		
	#if type(dataVersion) == unicode:
	#	dataVersion = dataVersion.encode('utf-8')
		
	#if type(tumorTypeSummary) == unicode:
	#	tumorTypeSummary = tumorTypeSummary.encode('utf-8')
	
	opened_output.writerow([str(geneSummary), str(highestSensitiveLevel), str(dataVersion), str(tumorTypeSummary), str(highestPrognosticImplicationLevel), str(oncogenic), str(alleleExist), str(citations), str(knownEffect), str(proteinStart), str(proteinEnd), str(tumorType), str(hugoSymbol),
	str(alteration), str(entrezGeneId), str(hgvs), str(consequence), str(referenceGenome), str(hotspot), str(diagnosticImplications), str(prognosticSummary), str(lastUpdate), str(prognosticImplications), str(treatments), str(geneExist), str(variantExist),
	str(highestDiagnosticImplicationLevel), str(diagnosticSummary), str(otherSignificantSensitiveLevels), str(otherSignificantResistanceLevels), str(variantSummary), str(vus), str(highestResistanceLevel)])
	
