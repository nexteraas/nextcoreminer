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


def write_output_peptide_aa(peptide_aa_dict, rank_df, header_list, outfile_name = 'data_cleaning_peptide_output.txt'):
	""" write the output file as an txt file for possible further manual processesing """
	
	with open(outfile_name, 'w') as of:
		## write the header preceeded by a # in order for easy removal
		for elmt in header_list:
			of.write(''.join(("#", elmt)))
			of.write('\n')

		## Write column names out
		of.write('AA_Peptide\tRank\tFrequency\tTotal_Count\tSource_counts\tSource_seq\tSource_list\n')
		## write the pandas df out
		for index, row in rank_df.iterrows():
			peptide_AA = row['peptide']
			tmp_source_count = ""
			tmp_source_seq = ""
			tmp_source_list = ""
			total_count = 0
			rank = str(row['rank'])
			freq = str(row['freq'])
			## make sure that the files are written out
			for indx, name in enumerate(peptide_aa_dict[peptide_AA]['seq_source_list']):
				## if there are more than one seq_source_list the index will be true >0
				if indx:
					## add the total count of the all the source sequnces 
					total_count += len(peptide_aa_dict[peptide_AA]['seq_list'][indx])
					## add it as a list to the dict but as eperate lists
					tmp_source_count = ','.join((tmp_source_count, str(len(peptide_aa_dict[peptide_AA]['seq_list'][indx]))))
					tmp_source_seq =  ','.join((tmp_source_seq, name))
					tmp_source_list = ','.join((tmp_source_list, str(peptide_aa_dict[peptide_AA]['seq_list'][indx])))
				else:
					total_count = len(peptide_aa_dict[peptide_AA]['seq_list'][indx])
					tmp_source_count =  str(len(peptide_aa_dict[peptide_AA]['seq_list'][indx]))
					tmp_source_seq = name
					tmp_source_list = str(peptide_aa_dict[peptide_AA]['seq_list'][indx])
			
			## Make the element for writing out the data
			temp_print = '\t'.join((peptide_AA, rank, freq, str(total_count), tmp_source_count, tmp_source_seq, tmp_source_list))
			
			of.write(temp_print)
			of.write('\n')

def write_output_peptide_xli(peptide_aa_dict, rank_df, header_list, outfile_name = 'data_cleaning_peptide_xli.txt'):
        """ write the output file as an txt file for possible further manual processesing """

        with open(outfile_name, 'w') as of:
                ## write the header preceeded by a # in order for easy removal
                for elmt in header_list:
                        of.write(''.join(("#", elmt)))
                        of.write('\n')

                ## Write column names out
                of.write('AA_Peptide\tRank\tFrequency\tTotal_Count\n')
                ## write the pandas df out
                for index, row in rank_df.iterrows():
                        peptide_AA = row['peptide']
                        tmp_source_count = ""
                        tmp_source_seq = ""
                        tmp_source_list = ""
                        total_count = 0
                        rank = str(row['rank'])
                        freq = str(row['freq'])
                        ## make sure that the files are written out
                        for indx, name in enumerate(peptide_aa_dict[peptide_AA]['seq_source_list']):
                                ## if there are more than one seq_source_list the index will be true >0
                                if indx:
                                        ## add the total count of the all the source sequnces
                                        total_count += len(peptide_aa_dict[peptide_AA]['seq_list'][indx])
                                        ## add it as a list to the dict but as eperate lists
                                        tmp_source_count = ','.join((tmp_source_count, str(len(peptide_aa_dict[peptide_AA]['seq_list'][indx]))))
                                        tmp_source_seq =  ','.join((tmp_source_seq, name))
                                        tmp_source_list = ','.join((tmp_source_list, str(peptide_aa_dict[peptide_AA]['seq_list'][indx])))
                                else:
                                        total_count = len(peptide_aa_dict[peptide_AA]['seq_list'][indx])
                                        tmp_source_count =  str(len(peptide_aa_dict[peptide_AA]['seq_list'][indx]))
                                        tmp_source_seq = name
                                        tmp_source_list = str(peptide_aa_dict[peptide_AA]['seq_list'][indx])

                        ## Make the element for writing out the data
                        temp_print = '\t'.join((peptide_AA, rank, freq, str(total_count))) #, tmp_source_count, tmp_source_seq, tmp_source_list))

                        of.write(temp_print)
                        of.write('\n')

def write_json_output(peptide_aa_dict, rank_df, header_list, outfile_name = 'data_cleaning_peptide_output.json'):
	""" write out the data in json format in order for further processing of the data"""
	
	peptide_out_dict = {}
	## Add the header as 'header' in the pre-json dict
	peptide_out_dict['header'] = header_list
	
	## iterate through the data frame
	for index, row in rank_df.iterrows():

			peptide_AA = row['peptide']

			peptide_out_dict[peptide_AA] = {}

			for indx, name in enumerate(peptide_aa_dict[peptide_AA]['seq_source_list']):
				## if there are more than one seq_source_list the index will be true > 0
				if indx:
					## add the total count of the all the source sequnces 
					peptide_out_dict[peptide_AA]['total_count'] += len(peptide_aa_dict[peptide_AA]['seq_list'][indx])
					peptide_out_dict[peptide_AA]['source_count'].append(str(len(peptide_aa_dict[peptide_AA]['seq_list'][indx])))
					peptide_out_dict[peptide_AA]['source_seq'].append(str(peptide_aa_dict[peptide_AA]['seq_source_list'][indx]))

				else:

					peptide_out_dict[peptide_AA]['rank'] = str(row['rank'])
					peptide_out_dict[peptide_AA]['freq'] = str(row['freq'])
					
					peptide_out_dict[peptide_AA]['total_count']  = len(peptide_aa_dict[peptide_AA]['seq_list'][indx])
					peptide_out_dict[peptide_AA]['source_count'] =  [str(len(peptide_aa_dict[peptide_AA]['seq_list'][indx]))]
					peptide_out_dict[peptide_AA]['source_seq'] = [str(peptide_aa_dict[peptide_AA]['seq_source_list'][indx])]
	
	
	## dump the dict to a json
	peptide_out_json = json.dumps(peptide_out_dict)

	## write out the json 
	with open(outfile_name, 'w') as of:
		of.write(peptide_out_json)

def data_list_handler_json(files_holder_list, n_process):
	## LOAD DATA
	### Do this in paralell in order to speed this up
	pepide_dict_list = []
	with concurrent.futures.ProcessPoolExecutor(max_workers = n_process) as executor:
		futures = []

		for file in files_holder_list:

			future = executor.submit(read_file_json, file)

			futures.append(future)

		for future in concurrent.futures.as_completed(futures):
	
			pepide_dict = future.result()

			pepide_dict_list.append(pepide_dict)

	return pepide_dict_list

def read_file_json(file_in):
	# open the file	
	with open(file_in, 'r') as f:		
		return json.load(f)

def data_for_json(peptide_dict, rounds_list, desitnation_folder):
	""" collection of the data needed to print out the jsonp file, for the data table visualisation """
	## take the rounds and sort them alpahnumerically for a nicer looking data table

	sorted_rounds_list = sort_on_source(rounds_list)
	
	
	## the informaiton used in the data tabel
	print_list = ['rank', 'freq', 'count', 'meta']

	pre_pd_dict = {}

	max_value_dict = {}

	R0_dict = {}

	for peptide in peptide_dict:

		## Since R0 is a non speceifc panning, ignore all peptides exclusivily found in this round 
		## This is done to minimize the amount of data that needs to be loaded
		## These are not discarted but saved in a seprate file 
		not_only_R0_flag = check_rounds(peptide_dict[peptide].keys())
		
		## if found in more than R0
		if not_only_R0_flag > 0:
			
			## For that peptide initate the dict
			pre_pd_dict[peptide] = {}
			
			## sorted_rounds_list is used to get the rigth column formatting out
			for _round in sorted_rounds_list: 
				##print_list is used since only a subset of the columns provided is used
				for attr in print_list:

					## Make a column name for the run and attriburte for easy distinction
					round_attr = '_'.join((_round, attr))
					## initaite the list instance
					if _round in peptide_dict[peptide]:
						
						## get the max count for the round in order to force the lowest scoring to the bottom
						## this is done for formatting, and sorting purposes in the data table 
						if round_attr in max_value_dict  and re.search(r'meta',round_attr) == None:
						
							if re.search(r'count', round_attr):
								if max_value_dict[round_attr] > peptide_dict[peptide][_round][attr]:
									max_value_dict[round_attr] = peptide_dict[peptide][_round][attr]
							else:
								if max_value_dict[round_attr] < peptide_dict[peptide][_round][attr]:
									max_value_dict[round_attr] = peptide_dict[peptide][_round][attr]
						else:
							max_value_dict[round_attr] = peptide_dict[peptide][_round][attr]
						
						## if the round is found add the attribute
						pre_pd_dict[peptide][round_attr] = peptide_dict[peptide][_round][attr]

					else:
						## if the round is not found add NA for easy subsitution later
						pre_pd_dict[peptide][round_attr] = "NA"
		else:
			## Add the ignored peptide to the R0 dict
			R0_dict[peptide] = peptide_dict[peptide]

	
	## generate the file paths for the data output
	data_output_destination = os.path.join(desitnation_folder, "data_table.json")
	R0_only_output_destination = os.path.join(desitnation_folder, "R0_data_table.json")

	## write the data output files to jsonp format
	write_out_json(pre_pd_dict, max_value_dict, data_output_destination)
	write_out_json(R0_dict, max_value_dict, R0_only_output_destination)

	## THE RETURN STATEMENT CAN BE DELETED AT A LATER TIME
	# return pre_pd_dict, R0_dict, max_value_dict

def sort_on_source(rounds_list):
	"""
		Sort the data based on the rounds (R0-R#), but collected by data source
	"""
	## HACK to get R0 to be the frist instance
	sorted_rounds = sorted(rounds_list)
	
	sorted_rounds_list = []
	round_id_indx = {}
	
	for indx, round_name in enumerate(sorted_rounds):
		if re.match(r'R0',round_name):
			sorted_rounds_list.append(str(round_name))
		else:
			## If not R0 check to see if it has been encounted before
			round_id = re.sub(r'R[1-9]_','',round_name)

			if round_id in round_id_indx:
				round_id_indx[round_id].append(indx)
			else:
				round_id_indx[round_id] = [indx]

	# Append all the rounds together based on the data source
	for _round in round_id_indx:
		for indx in round_id_indx[_round]:
			sorted_rounds_list.append(sorted_rounds[indx])

	return sorted_rounds_list

def check_rounds(rounds_list):
	""" Counts the round that are not R0"""
	round_count = 0
	for _round in rounds_list:

		if re.search(r'^R0', _round) is None:		
			round_count += 1			

	return round_count

def write_out_json(data_in, max_value_dict, desitnation):

	data_in_outformat = data_structure(data_in, max_value_dict)
	
	with open(desitnation,'w') as of:
		data_file_out_json = json.dumps(data_in_outformat)	
		of.write(data_file_out_json)

def data_structure(data_in, max_value_dict):
	""" 
		Save all the data in the correct format for data tables 
		Subsitute all NA with the appropriate max value provided by the max_value_dict
		This is done in a list of lists that can be paresed by the javascript application

		Since all peptides contain either a value of NA for all attributes the NA's can easly be substituted
		This also enables a 1-time runthough of the dict
	"""
	
	data_list = []
	
	for peptide in data_in:
		
		data_dict = {}
		## save all the unique peptides
		data_dict['peptide'] = peptide


		for attr in data_in[peptide]:
			## if the attribute is NA take the approriate action
			if data_in[peptide][attr] == 'NA':
				## all actions will force the given peptide to the bottom if sorted on the given attr. column
				if re.search(r"_freq",attr):
					data_dict[attr] = 0
				elif re.search(r"_count",attr):
					if attr in max_value_dict:
						data_dict[attr] = int(max_value_dict[attr])*-1
					else:
						data_dict[attr] = -1

				elif re.search(r'_meta', attr):
					data_dict[attr] = [['NA'],['NA']]
				else:
					if attr in max_value_dict:
						fltAttr = int(max_value_dict[attr])
						data_dict[attr] = fltAttr*100
					else: 
						data_dict[attr] = 100000

			## else add the infromation to the data structure
			else:
				data_dict[attr] = data_in[peptide][attr]

		#append the data to the list to make list olists
		data_list.append(data_dict)

	return data_list

def generate_workfolder(NNK_name, parent_folder):
	""" 
		Generate a output folder for each experiment run
		This is done in order to ensure that data is not overwritten at any point
		The naming is a combination of the time of day and a a set of three radom numbers and three random letters
		Since there is a check fo an excitsing name ther should not be any cases of overwritten files
		evne if the files are run at the exact same time
	"""
	## time of the run
	today = datetime.datetime.now()
	date_today = today.strftime("%Y%m%d")
	## generaton to the three letter, three number code for each run
	rand_name = generate_random_name()
	rand_name_out = today.strftime("%H%M%S") +"_"+ rand_name
	
	run_dir = os.path.join(parent_folder,"Data", NNK_name, date_today, rand_name_out)
	web_dir = os.path.join(NNK_name, date_today,rand_name_out)
	## check if the run/folder excists
	if not os.path.exists(run_dir):
		## if not create all the folders needed to generate the full path
		os.makedirs(run_dir)
		## otherwise generate a new name 
	else:
		rand_name = generate_random_name()
		rand_name_out = today.strftime("%H%M%S") +"_"+ rand_name
		run_dir = os.path.join(parent_folder,"Data", NNK_name, date_today, rand_name_out)	
		web_dir = os.path.join( NNK_name, date_today,rand_name_out)
		os.makedirs(run_dir)

	return(run_dir,web_dir)

def generate_random_name():
	""" generation of the random three letter, three number tage"""
	rand_char = ''.join([random.choice(string.ascii_letters) for i in range(3)])
	rand_num = ''.join([str(random.randint(0, 9)) for i in range(3)])
	name = ''.join([rand_char, rand_num])
	return name
