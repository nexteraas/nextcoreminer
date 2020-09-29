#!/usr/bin/python3

import argparse
import os
import re
import sys
from pathlib import Path

import concurrent.futures

import datetime 
import numpy as np
import pandas as pd
import random
import string

import json

def dna_pwm_regex_gen(regex_string, split= 1):
	"""	
	Construction of the Postional Weigted Matrix (PWM), needed for sequence selection, based on the library regex provided from the 
	mandatory user input args.dnaregexstring (-d). How this parameter is processed is based on the Galaxy input.
	And migth need to be changed when we design the input from the library in Galaxy
	"""
	
	regex_string_ref_dict = dict()

	if split:
		# migth need to be changed when we design the input from the library in Galaxy, since this  mi
		regex_string_edit = regex_string.split('_')	

	else:
		regex_string_edit = regex_string
	
	# Hand made reference dict for translation of the regex
	regex_ref_dict = {
	'A':'A',
	'T':'T',
	'C':'C',
	'G':'G',
	'A':'A',
	'M':'CA',
	'Y':'TC',
	'K':'TG',
	'W':'AT',
	'R':'AG',
	'B':'TCG',
	'D':'ATG',
	'H':'ATC',
	'V':'ACG',
	'N':'ATCG'
	}
    
	n_p = 0
	prob_base_0 = 1
	
	# Constuct a PWM in order to do a sum of any match made
	for codon in regex_string_edit:
		for nuc in codon:
			regex_string_ref_dict[n_p] = dict()
			regex_string_ref_dict[n_p]['A'] = prob_base_0
			regex_string_ref_dict[n_p]['T'] = prob_base_0
			regex_string_ref_dict[n_p]['C'] = prob_base_0
			regex_string_ref_dict[n_p]['G'] = prob_base_0
			regex_string_ref_dict[n_p]['N'] = prob_base_0

			for ref_nuc in regex_ref_dict[nuc]:
				regex_string_ref_dict[n_p][ref_nuc] = 0
			n_p += 1

	return regex_string_ref_dict

def add_to_dict(dict_peptide, peptide, sample_name):
	# Add any sample to a like peptide entry and update the source count of that entry 
	# If the entry dose not exist make it

	if peptide not in dict_peptide:
		dict_peptide[peptide] = {}
		dict_peptide[peptide]['sample_name'] = sample_name
	else:
		dict_peptide[peptide]['sample_name'] + sample_name
	
	dict_peptide[peptide]['sample_count'] = len(dict_peptide[peptide]['sample_name'])

	return dict_peptide
def extract_barcode(filename):

	filename_array = filename.rstrip('\r\n').split('_')

	if re.match(r"^([ATCGN])+$", filename_array[3]):
		return filename_array[3]
	else:
		for elm in filename_array:
			if re.match(r"^([ATCGN])+$",elm):
				return elm

def read_to_dict(filename, barcode):

	"""
	Read through the sequnece file, and extract the sequences for which the bar-code provided can be found.
	If no barcode is found, read through the filename and try to extract the bar-code
	Discard all the reads that do no contain correct bar-code, in the sequence identification line, and in the sequence
	cut the sequence at the end of the bar-code, and get the peptide
	"""	

	## For summary stats 
	seq_count = 0

	#The output dict construction
	dict_nuc_seq = {}
	dict_barcode = {}

	flag = 0 

	with open(filename, 'r') as f:

		for line in f:
			
			_line = line.rstrip('\r\n')

			if flag:
				##When id line i found, the next line is the sequence, set the flag to 0, as there is only one sequence line
				flag = 0

				##Make sure that it is a DNA sequence read
				if re.match(r"^([ATCGN])+", _line):
					seq_count += 1
					
					##Find the barcode, if not found return -1 
					indx_barcode = _line.find(_barcode)
					##if barcode not found skip the line
					if indx_barcode != -1:	
						
						##Extract the infomation after the end of the barcode
						_line_ = _line[indx_barcode+len(_barcode):]
						
						##if the extraction dose not retun an empty sting move on
						### Can be set to an interger in length of expected peptide
						if _line_ != '':

							## Look for the peptide, and make sure that if found that the barcode are the same
							if _line_ in dict_nuc_seq and dict_nuc_seq[_line_]['barcode'] == _barcode:
								##if true add to the dict
								dict_nuc_seq[_line_]['seq_list'].append(name)
							else:
								##else construct a new enty in the dict
								dict_nuc_seq[_line_] = {}
								dict_nuc_seq[_line_]['seq_list'] = [name]
								dict_nuc_seq[_line_]['barcode'] = _barcode

			#if the identification line, denoted by @ in the start if found, set the name to line and set acceptance flag to 1
			if re.search(r'^\@\w+:', _line.rstrip('\r\n')):
				##split on ':' in order to extract the barcode
				name_line_array  = _line.split(':')
				##The @header must end with the barcode in order for extraction
				if re.match(r"^([ATCGN])+$",name_line_array[len(name_line_array)-1]):
					_barcode = name_line_array[len(name_line_array)-1]
				else:
					_barcode = barcode
				
				name = _line
				flag = 1
				## if the barcode is known add a count to the barode
				if _barcode in dict_barcode:
					dict_barcode[_barcode] += 1
				## else add the entry
				else:
					dict_barcode[_barcode] = 1

	## Extract the most abundant barcode and use this	
	barcode_to_use = select_barcode(dict_barcode)

	### Select the sequences that orignate from the most abundant barcode
	dict_nuc_seq = select_seq_based_on_barcode(dict_nuc_seq, barcode_to_use)
	
	return dict_nuc_seq, barcode_to_use, seq_count, dict_barcode

def select_barcode(barcode_dict):
	"""Select the most abundant barcode by count"""
	
	selection_threshold = 0
	barcode_to_use = ''
	for barcode in barcode_dict:
		if selection_threshold < barcode_dict[barcode]:
			selection_threshold = barcode_dict[barcode]
			barcode_to_use = barcode

	return barcode_to_use

def select_seq_based_on_barcode(dict_in, barcode):
	"""Select the sequneces that originate from the most abundant barcode, return these"""
	dict_out = {}
	
	for key in dict_in:
		if barcode == dict_in[key]['barcode']:
			dict_out[key] = dict()
			for info_key in dict_in[key]:
				dict_out[key][info_key] = dict_in[key][info_key]

	return dict_out

def get_get_substring(seq, size, pwm_dict, max_mismatch):
	""" Match, count mismathces and retrun a substing of the sequncen based on the PWM provided"""
	
	## Return obejcts
	suqbstrings_match = ''
	best_count = None

	##iterate over all positions of the sequence
	for i in range(len(seq)):
		## If only a parital peptide can be obtained, break
		if i + size > len(seq):
			break
		else:
			## Subset the sequence to be matched
			subset_seq = seq[i:(i+size)]
			
			## get the mismatch count when comparing the sequnece to the PWM
			count = compare_to_pwm(subset_seq, pwm_dict, max_mismatch)

			if count is not None:
				## best count will only be None if it is the first comparison
				if best_count is not None:
					## if the new count is lower
					if best_count > count:
							### not used at the moment, can be re-impelmented in order the allow for ties
							## else frist-best is used
							if best_count == count:
								print('tie')
								suqbstrings_match = subset_seq
							#if ther is a lower mismatch count, then set the new substring
							else:
								suqbstrings_match = subset_seq
							## update the best_count for the match
							best_count = count

				## For the first comparison set the best_count and and the subset_match
				else:
					best_count = count
					suqbstrings_match = subset_seq
	## if the substing is not empty then retrun
	if suqbstrings_match:
		return suqbstrings_match, best_count
	else:
		return None, None

def add_to_dict_out_2(dict_peptide, peptide, dict_nuc_seq, sequence):
	
	if peptide not in dict_peptide:
		dict_peptide[peptide] = dict()
		dict_peptide[peptide]['seq_list'] = dict_nuc_seq[sequence]['seq_list']
	else:
		dict_peptide[peptide]['seq_list'] +  dict_nuc_seq[sequence]['seq_list']
		
	return dict_peptide

def compare_seq_to_pwm(peptide_dict, pwm_dict, max_mismatch):
	""" Comaper the sequnece a given PWM, get the best matching substing and the mismatch count"""
	dna_pepetide_dict = {}
	## get each reading frame/sequence in the sequnce analysis run
	for peptide_rf in peptide_dict:
		
		## pass the full sequnece and compare to the PWM, and get the peptide reading frame back
		peptide_dna, count = get_get_substring(peptide_rf, len(pwm_dict), pwm_dict, max_mismatch)
		
		## if the return is not None add it tho the reading frame dict
		if peptide_dna:
			## if it is already found then add it to the exsiting dict
			if peptide_dna in dna_pepetide_dict:
				dna_pepetide_dict[peptide_dna]['mismatch_count'].append(count)
				dna_pepetide_dict[peptide_dna]['seq_list'].extend(peptide_dict[peptide_rf]['seq_list'])
			## else make the enty
			else:
				dna_pepetide_dict[peptide_dna] = {}
				dna_pepetide_dict[peptide_dna]['mismatch_count'] = [count]
				dna_pepetide_dict[peptide_dna]['seq_list'] = peptide_dict[peptide_rf]['seq_list']


	return dna_pepetide_dict

def compare_to_pwm(seq, pwm_dict, max_mismatch):	
	""" Comaper the sequnece a given PWM, and is the mismatchs do not excede the max-mismatch count return the count"""
	count = 0
	## loop over location and nuc
	for i, nuc in enumerate(seq):
		## check if it is the PWM
		if pwm_dict[i][nuc]:
			count += pwm_dict[i][nuc]
		if max_mismatch < count:
			return None
	return count


def search_by_known_flank(seq_dict, upstream, downstream, pwm_siz,  max_mismatch):
	""" Compare the sequende to the flanks passed, and extract the sequence between these """
	print(upstream, downstream)
	discarted_seq = {}
	## Generate a PWM for each flanks
	pwm_upstream = dna_pwm_regex_gen(upstream, 0)
	pwm_downstream = dna_pwm_regex_gen(downstream, 0)

	peptide_found_dict = dict()
	
	### SUMMARY STATS VARIABLES
	downstream_found = 0
	upstream_found = 0
	all_seq  = 0 
	
	upstream_manual = 0
	downstream_manual = 0

	rf_summary_dict = dict()

	## Loop over all the unique sequences extracred from the data soruce
	for sequence in seq_dict:

		## Compare the flank the the sequnces, return the index of the best fit
		upstream_index, upstream_count = compare_flank_pwm(sequence, pwm_upstream, max_mismatch)
		## continue the only it the upstream indes is found
		if upstream_index is not None:
			## subset the sequnce, removing the uptream flank
			upstream_flank_seq = sequence[upstream_index + len(upstream): ]
			## Compare the flank the the sequnces, return the index of the best fit
			downstream_index, downstream_count = compare_flank_pwm(upstream_flank_seq , pwm_downstream, max_mismatch)
			## continue the only it the downtream indes is found
			if downstream_index is not None:

				## set the base of the downstream_index in order to get the correct out sequence
				downstream_index += upstream_index + len(upstream)

				## extract the peptide reading sframe
				peptide_rf  = sequence[upstream_index + len(pwm_upstream) : downstream_index]

				## For each reading frame length add the lenth to a sum stat dict
				if str(len(peptide_rf)) in rf_summary_dict:
					rf_summary_dict[str(len(peptide_rf))] += 1
				else:
					rf_summary_dict[str(len(peptide_rf))] = 1

				## Only allow reading frames that are dicisable by 3 and of the minimum length of the PWM
				## The peptide generation in the source data dictates that all generated peptides are codon based and will therefore be divisiable by 3
				if len(peptide_rf) % 3 == 0 and len(peptide_rf) / pwm_siz >= 1:

					if re.match(r'[N]',peptide_rf):
						print(peptide_rf)
						discarted_seq[peptide_rf] = {}
						discarted_seq[peptide_rf][sequence] = len(seq_dict[sequence]['seq_list'])

					else:
						## if the reading frame is an entry in the dict add the sequence
						if peptide_rf in peptide_found_dict:
				
							peptide_found_dict[peptide_rf]['seq_list'].extend(seq_dict[sequence]['seq_list'])
						## else contruct the enty
						else: 
							
							peptide_found_dict[peptide_rf] = dict()
							peptide_found_dict[peptide_rf]['seq_list'] = seq_dict[sequence]['seq_list']

	return peptide_found_dict, rf_summary_dict, discarted_seq


def compare_flank_pwm(seq, pwm_dict, max_mismatch):
	""" Compare the flanking PWM to the provided sequnce get the index of the best macth, return these"""	
	min_count = 1000
	min_start = None

	## for each position in the sequence comapre pos -> len(pwm_dict)
	for pos in range(len(seq)):
		
		#assume that there is a full match	
		pwm_start = 0
		
		min_count_inner = 1000
		min_start_inner = None
		
		## in order to allow a not full match add a penalty to the count
		## As long as this do not exceed the max_mismatch continue the comparion
		while pwm_start <= max_mismatch:
			
			## set the internal conuter to 0 at the start of the comparison
			count = 0
			## compare the length of the PWM 
			for i in range(pwm_start, len(pwm_dict)):
				
				## if the start is shifted add the apporopiate penalty
				count += pwm_start
				
				index = pos + i - pwm_start
				
				if index < len(seq):
					count += pwm_dict[i][seq[index]]		
				
				else:
					count += abs(index + len(pwm_dict) - i -len(seq))
					break

			## Update the counter for the the shift
			if min_count_inner > count and count <= max_mismatch:
				
				min_count_inner = count
				min_start_inner = pos

			pwm_start += 1

		if min_count > min_count_inner and min_count_inner <= max_mismatch:
			min_count = min_count_inner
			min_start = min_start_inner

	if min_count <= max_mismatch:
		return min_start, min_count
	else: 
		return None, None 

def print_reading_frame_stats(rf_summary_dict):
	""" Print function for the reading frame lenght dict"""
	
	print('reading_frame_length\tsequneces count')

	keys_order = sorted([int(i) for i in rf_summary_dict])
	for i in keys_order:
		print('\t'.join((str(i), str(rf_summary_dict[str(i)]))))
	print('\n')

def print_discarted_preptide_rf(discarted_seq):

	for reading_frame in discarted_seq:
		print(discarted_seq)
		for sequnece in discarted_seq[reading_frame]:
			print("\t".join([reading_frame, sequnece, discarted_seq[reading_frame][sequnece]]))
			print("\n")

def do_paralell_translation(peptide_dna_seq_list, translation_lib_reference, n_process=6):	
	""" In order to speed up the procsses if translation, do it in parallel
		This if there is not a translation libray provided then use the universel 
		one providede, for the translation
	"""
	
	peptide_dict_out = {}
	
	### if there is no translaton lib provided then use this
	if translation_lib_reference == None:
	    print('Codon Table used: Univesal')
	    codon_table = { 
	        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
	        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
	        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
	        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
	        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
	        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
	        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
	        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
	        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
	        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
	        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
	        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
	        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
	        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
	        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
	        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
	    } 
	else:
		## else get the table given
		codon_table = get_translation_codon_table(translation_lib_reference)


	## multiprocessing	
	with concurrent.futures.ProcessPoolExecutor(max_workers = n_process) as executor:
	
		futures = {}
		
		for _peptide in peptide_dna_seq_list:

			future = executor.submit(translate_dna_to_protein, _peptide, codon_table)
			
			futures[future] = _peptide

		## uncovelute the runs
		for future in concurrent.futures.as_completed(futures):
			
			peptide_aa = future.result()

			ref_peptide = futures[future]
			
			## if there is a return translated peptide then add it as an new enty or to an entry
			if peptide_aa:
				## if translated peptide is already an entry the add the new source
				if peptide_aa in peptide_dict_out:

					peptide_dict_out[peptide_aa]['seq_list'].append(peptide_dna_seq_list[ref_peptide]['seq_list'])
					peptide_dict_out[peptide_aa]['seq_source_list'].append(ref_peptide)
				## else create it as an entry
				else:
					peptide_dict_out[peptide_aa] = {}

					peptide_dict_out[peptide_aa]['seq_list'] = [peptide_dna_seq_list[ref_peptide]['seq_list']]
					peptide_dict_out[peptide_aa]['seq_source_list'] = [ref_peptide]

	return peptide_dict_out

def get_translation_codon_table(condon_table_ref_name):
	""" THE CODON TABLES MUST BE LOCATED IN ~/codon_ref_lib/* otherwise the progam will not work"""
	json_file_name = '/home/Nextera/galaxy/tools/testing/codon_ref_lib/' + str(condon_table_ref_name) + '.json'
	try:
		with open(json_file_name, 'r') as f:
			codon_ref_table = json.load(f)
			return codon_ref_table
	except ValueError:
		print("The codon reference table could not be read")
		print("Control if the files are located in the '~/codon_ref_lib/' folder and input table name")

def translate_dna_to_protein(peptide, translation_matrix = None):
	""" Translate the DNA peptide into the the corrosponding AA translated peptide """
	codon_table = translation_matrix

	peptide_aa = ''
	try:
		# ensure that the assumption of the peptide generation in the source data dictates that all generated peptides are codon based and will therefore be divisiable by 3
		if len(peptide)%3 == 0: 
			## take the sequence the in codon chunks
			for i in range(0, len(peptide), 3): 
				## translate the codon
				codon = peptide[i:i + 3] 
				## ensure that the codon is DNA
				if re.match(r'[A-Za-z]', codon_table[codon]) is not None:
					## Add the translated condon to the sequence
					peptide_aa += codon_table[codon] 
				else:
					return None
		## retun the translated peptide
		return peptide_aa
	except ValueError:
		print("Invalid peptide length\nPeptide not divisible by 3")
  
def print_output_peptide_aa(peptide_aa_dict):
	"""Print function for the translated peptide dict"""
	print('AA_Peptide\tTotal_Count\tSource_counts\tSource_seq\tSource_list\n')
	for peptide_AA in peptide_aa_dict:
		tmp_source_count = ""
		tmp_source_seq = ""
		tmp_source_list = ""
		total_count = 0

		for indx, name in enumerate(peptide_aa_dict[peptide_AA]['seq_source_list']):
	
			if indx:
				total_count += len(peptide_aa_dict[peptide_AA]['seq_list'][indx])
				tmp_source_count = ','.join((tmp_source_count, str(len(peptide_aa_dict[peptide_AA]['seq_list'][indx]))))
				tmp_source_seq =  ','.join((tmp_source_seq, name))
				tmp_source_list = ','.join((tmp_source_list, str(peptide_aa_dict[peptide_AA]['seq_list'][indx])))
			
			else:
				total_count = len(peptide_aa_dict[peptide_AA]['seq_list'][indx])
				tmp_source_count =  str(len(peptide_aa_dict[peptide_AA]['seq_list'][indx]))
				tmp_source_seq = name
				tmp_source_list = str(peptide_aa_dict[peptide_AA]['seq_list'][indx])
		
		temp_print = '\t'.join((peptide_AA, str(total_count), tmp_source_count, tmp_source_seq, tmp_source_list))
	
		print(temp_print)

def rank_peptides(peptide_aa_dict):
	""" Add a rank to the translated peptide
		Store the data in a pandas data frame, in order to rank and sort the data
		Ranking is done by the following rules:
			Highest count gets best rank [1]
			If count is a tie, the lowest number of unique source sequences gets the best rank
			If count is a tie, and the number of unique source sequnecese are tie give same rank
		retun the pandas data frame
	"""

	## Define column names of the pandas dataframe
	cols = ['peptide', 'count', 'source', 'rank']
	data_list = []
	total_count = 0
	## loop over all the translated peptide
	for peptide_AA in peptide_aa_dict:
		_count = 0
		_source_numb = 0
		## for each list in the list of source sequences for the translated sequence 
		for indx, name in enumerate(peptide_aa_dict[peptide_AA]['seq_source_list']):
			## Get the count for each list, and add that to the full count for that translate
			_count += len(peptide_aa_dict[peptide_AA]['seq_list'][indx])
			total_count += len(peptide_aa_dict[peptide_AA]['seq_list'][indx])
		
		## Get the number of unique sequneces that lead to the construction of the translated peptide
		_source_numb = len(peptide_aa_dict[peptide_AA]['seq_source_list'])
		## Append the this list to the the for translation into pd
		data_list.append([peptide_AA, _count, _source_numb, 0])

	## make the data into a data frame
	df_ = pd.DataFrame(data_list, columns=cols)
	
	## sort values by count in in and decending order and the numver of sources in an acending order
	df_ = df_.sort_values(by=['count', 'source'], ascending=(False, True))
	## envoke the new sorting by droping the index
	df_ = df_.reset_index(drop=True)
	
	running_count = 0
	running_rank = 1
	running_source_count = 0
	## ASSING THE RANK BASED ON COUNT THEN SOURCES
	for index, row in df_.iterrows():
		#if the count is a lower than the count before add one to rank and assign the rank
		
		if running_count > row['count']:		
			
			running_rank += 1
			running_count = row['count']
			running_source_count = row['source']
			
			df_.at[index, 'rank'] = running_rank

		#if the count is equal check the source count
		elif running_count == row['count']:
			#if the sequnce count is higher then add one to the rank	
			
			if running_source_count < row['source']:
				running_rank += 1
				running_source_count = row['source']
				df_.at[index, 'rank'] = running_rank
			else:
				df_.at[index, 'rank'] = running_rank
		
		if running_rank == 1:
			df_.at[index, 'rank'] = running_rank
			running_count = row['count']
			running_source_count = row['source']
	
	## add the frequency 
	df_['freq'] = df_['count'].div(total_count).round(10)
	## sort by rank 
	df_ = df_.sort_values(by='rank', ascending=True)
	## envoke the new sorting by droping the index
	df_ = df_.reset_index(drop=True)

	return df_


def get_sequnce_count(in_dict, nested = False):
	""" counter for the sum stat dict """
	count = 0
	for i in in_dict:
		if nested == 1:
			for indx, name in enumerate(in_dict[i]['seq_list']):
				
				count += len(in_dict[i]['seq_list'][indx])
		else:
			count += len(in_dict[i]['seq_list'])

	return count

def generate_header(soruce_filename, barcode, flanks, NNK_lib_name, *args):
	""" 
		Generate the header for the output files
		The information given in the *args needs to have a key:value pair in order to be added
	"""
	## get runtime in order to compare to the further collection and comparison of the data
	now = datetime.datetime.now()
	run_date = now.strftime("%Y%m%d:%H:%M:%S")
	header = [soruce_filename]
	##Find the index of round description in the file name
	print(soruce_filename)
	## Add the round in 
	round_index_start = re.search(r'(R[0-9]+)_([A-Za-z0-9]+)', str(soruce_filename)).start()
	round_index_end =  re.search(r'(R[0-9]+)_([A-Za-z0-9]+)', str(soruce_filename)).end()
	## concat round name
	round_name = soruce_filename[round_index_start:round_index_end]

	header.append(round_name)
	header.append(barcode)
	## add flanks for backtrace
	flanks_out = ':::'.join((flanks[0],flanks[1]))
	header.append(flanks_out)
	## Append NNK Lib for backtrace
	header.append("NNK_lib_used::"+ str(NNK_lib_name))
	header.append('run_date:' + run_date)
	## join the header, and retrun it 
	if args and type(args) == 'dict':
		for elemt in args:
			header.append(':'.join((str(elemt), str(args[elemt]))))
	return header	
		
def collect_peptide_dict_from_json(pepetide_dict_list):
	""" Collect the peptide infomation from the peptide dict list into a compairabel format for all the parssed studies"""
	
	collected_dict = {}
	rounds_present = []
	header_info = []
	

	## for each experiment in the list extract he relevant information 
	for peptide_dict in pepetide_dict_list:
		## this position is given in the data_clean.py script
		_round = get_round(peptide_dict['header'])
		## for each peptide in the given experiment
		for peptide in peptide_dict:	
			if peptide != 'header':
				
				if _round not in rounds_present:
					rounds_present.append(_round)
				## if peptide enty has already been made
				if peptide in collected_dict:
					## add add the new entry into the allready initated dict, and add the appropriate NA if any
					collected_dict = add_peptide_info_json(peptide,_round, peptide_dict, collected_dict)
					
				## otherwise initiate the entry
				else:
					collected_dict[peptide] = {}
					## add add the new entry into the allready initated dict, and add the appropriate NA if any
					collected_dict = add_peptide_info_json(peptide, _round, peptide_dict, collected_dict)
			## if header extract the information needed to make the appropriate folder structure
			else:
				NNK_name = get_header_info(peptide_dict['header'])
				header_info.append(NNK_name)		
	## Write out if the comapriason is made between different NNK libaries
	if len(set(header_info)) != 1:
		print("A comparison of differnet NNK_libraies are be made")
		print("Saved under loaction: " + header_info[0] + " folder loaction")
		print("NNK libraries found are:")
		for NNK in header_info:
			print(NNK)

	return collected_dict, rounds_present, header_info[0]

def get_round(header_array):

	for n, head_info in enumerate(header_array):
		if re.match(r'fastq', head_info) is None:
			if n != 0 and re.match(r'^R[0-9]_', head_info):
				return head_info

def get_header_info(header):
	""" extract the needed inforamtion from the header to be used in the the folder construction """
	NNK_tag_full = None
	for header_elm in header:
		if re.match(r'NNK_lib_used', header_elm):
			NNK_tag_full = header_elm

	NNK_tag_list = NNK_tag_full.split('::')

	NNK_tag = NNK_tag_list[len(NNK_tag_list)-1]

	return NNK_tag

def add_peptide_info_json(peptide, _round, peptide_dict, output_dict):
	""" add the data from the json file to the dict for the peptide entry """
	output_dict[peptide][_round] = {}
	output_dict[peptide][_round]['rank'] = peptide_dict[peptide]['rank']
	output_dict[peptide][_round]['freq'] = peptide_dict[peptide]['freq']
	output_dict[peptide][_round]['count'] = peptide_dict[peptide]['total_count']
	# Collect the DNA source and the counts from these together into two arrays
	output_dict[peptide][_round]['meta'] = [peptide_dict[peptide]['source_count'],peptide_dict[peptide]['source_seq']]
	output_dict[peptide][_round]['other'] = peptide_dict['header']

	return output_dict

def peptide_nice_fromat(peptide_dict, rounds_list):
	""""
		print out the peptide dict in a nice fromat for inspection
	"""
	sorted_rounds_list = sorted(rounds_list)
	print_list = ['rank', 'freq', 'count']
		
	header = ['peptide']
	for _round in sorted_rounds_list:
		for attr in print_list:
			header.append(_round+"_"+attr)

	print('\t'.join(header))


	for peptide in peptide_dict:
		
		peptide_information = peptide
		
		for _round in sorted_rounds_list:
			
			round_info = None
			
			for attr in print_list:
				if round_info is not None:
					if _round in peptide_dict[peptide]:
						
						round_info = '\t'.join((round_info, str(peptide_dict[peptide][_round][attr])))
					else:
						round_info = '\t'.join((round_info,'NA'))
				else:
					if _round in peptide_dict[peptide]:
						round_info = str(peptide_dict[peptide][_round][attr])
					else:
						round_info = 'NA'

			peptide_information = '\t'.join((peptide_information,round_info))
			
		print(peptide_information)


