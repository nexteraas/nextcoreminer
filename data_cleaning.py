#!/usr/bin/python3
import argparse
import os
import re
import sys
import pandas as pd
import numpy as np

import datetime 

import json

import concurrent.futures

import viz_lib as vl
import data_preprocessign_methods as dpm
import file_handler as fh


def parse_arguments():
	""" All arugments must be there in oder for the program to run"""

	parser = argparse.ArgumentParser()
	
	parser.add_argument('-s', '--seqfile', default=sys.stdin, required=True,
						help='Path to the input file.'
						'The format must be .fastq, since the read-through looks for the identifiers in the specific format'
						)
	parser.add_argument('-n', '--filename', default=sys.stdin, required=True)
	parser.add_argument('-d', '--dnaregexstring', default=sys.stdin, required=True,
						help='The library regular expression to the used in the PWM construction'
						)
	parser.add_argument('-r', '--NNK_lib_name', default=sys.stdin, required=True,
						help='The library name for the regular expression to the used in the PWM construction.'
						'Used in the header for the downstream folder generation at referencing'
						)
	parser.add_argument('-f', '--flanking_lib', default=sys.stdin, nargs='+',
						help='The flaking library used to identify the flaking region'
						)
	parser.add_argument('-g', '--flanking', default=sys.stdin,
						help='The number of maximum up and downstream flanking region mismatches'
						)
	parser.add_argument('-c', '--codontable', default=None,
						help='The name of the codon table to use when translating from DNA to AA'
						'Default: [Univesal]'
						)
	parser.add_argument('-o', '--output', default=sys.stdin, required=True, help="output file")
	return parser.parse_args()

def main():
	
	args = parse_arguments()

	if args:
	
		dict_data_retention_overview = {}
		dict_data_sumstats = {}
		
		#DATA GENERATE PWM
		pwm_dict = dpm.dna_pwm_regex_gen(args.dnaregexstring)
		#GET THE BARCODE
		barcode = dpm.extract_barcode(args.filename)
		#WRITE SEQUENCES TO DICT
		dict_seq, barcode, raw_seq_count, dict_barcode_count = dpm.read_to_dict(args.seqfile, barcode)
		## SUMSTATS ##
		dict_data_retention_overview['raw_count'] = raw_seq_count
		dict_data_retention_overview['barcode_found'] = dpm.get_sequnce_count(dict_seq)
		dict_data_sumstats['seq_count']	= len(dict_seq)
		#DATA READINGFRAME FOUND
		dict_peptide_rf, dict_readingframe, discarted_seq = dpm.search_by_known_flank(dict_seq, args.flanking_lib[0], args.flanking_lib[1], len(pwm_dict), int(args.flanking))

		## SUMSTATS ##
		dict_data_retention_overview['flanks_found'] = dpm.get_sequnce_count(dict_peptide_rf)
		dict_data_sumstats['seq_count_flanks']	= len(dict_peptide_rf)		
		###
		### IF any print the discarted sequnces
		if len(discarted_seq) > 0:
			dpm.print_discarted_preptide_rf(discarted_seq)
			
		dpm.print_reading_frame_stats(dict_readingframe)
		#DATA FIND BEST TRANSLATION FRAME


		dict_peptide_dna = dpm.compare_seq_to_pwm(dict_peptide_rf, pwm_dict, int(args.flanking))
		
		## SUMSTATS ##
		dict_data_retention_overview['regex_pwm_matched'] = dpm.get_sequnce_count(dict_peptide_dna)
		dict_data_sumstats['seq_count_peptide_dna']	= len(dict_peptide_dna)		
		##
		
		#DATA TRANSLATION 
		dict_peptide_aa = dpm.do_paralell_translation(dict_peptide_dna, args.codontable)
		
		## SUMSTATS ##
		dict_data_retention_overview['amino_acid'] = dpm.get_sequnce_count(dict_peptide_aa, 1)
		dict_data_sumstats['seq_count_peptide_aa']	= len(dict_peptide_aa)		
		##
		
		vl.data_overview_viz("graph_out", dict_data_sumstats, dict_data_retention_overview, dict_readingframe, dict_barcode_count)
		## SUMSTATS VIZ
		
		#DATA RANKING
		df_rank = dpm.rank_peptides(dict_peptide_aa)
		#DATA HEADER GENERATION
		if args.codontable != None:
			codontable = str(args.codontable).lower()
		else:
			codontable = 'Universal'

		list_header_for_run = dpm.generate_header(args.filename, barcode, args.flanking_lib, args.NNK_lib_name, {'codon_table':codontable})


		#DATA WRITE OUTPUT

		fh.write_json_output(dict_peptide_aa, df_rank, list_header_for_run)
		fh.write_output_peptide_aa(dict_peptide_aa, df_rank, list_header_for_run)
		fh.write_output_peptide_xli(dict_peptide_aa, df_rank, list_header_for_run)


if __name__ == "__main__":
	main()
