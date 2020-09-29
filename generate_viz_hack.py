#!/usr/bin/python3
import argparse
import os

def parse_arguments():
	""" All arugments must be there in oder for the program to run"""

	parser = argparse.ArgumentParser()
	
	#parser.add_argument('-h', '--html_template', default=sys.stdin, required=True, nargs='+',
	#					help='Path to the .html template'
	#					)
	parser.add_argument('-f', '--jsonp_File', default=4, required=True,
						help='Path to the .jsonp file'
						)
	#parser.add_argument('-j', '--javascript_template', default=4, required=True,
	#					help='Path to the javascript template'
	#					)
	return parser.parse_args()

def html_for_galaxy(html_template):

	with open('data_display.html', 'w') as of:
		with open(html_template,'r') as f:

			for line in f:
				of.write(line)

def write_javascritpfile(javascritp_template, jsonppos):

	with open('scritp.js', "w") as of:

		with open(javascritp_template, "r") as f:
			n = 1
			for line in f:

				if n == 65:
					of.write("url:"+jsonppos)
					of.write("\n")
				else:
					of.write(line)
def main():
	args = parse_arguments()

	if args:
		html_for_galaxy('/home/Nextera/galaxy/tools/testing/table_frontend_collection/data_display_template.html')
		write_javascritpfile('/home/Nextera/galaxy/tools/testing/table_frontend_collection/data_table_viz/js/script_template.js', args.jsonp_File)

if __name__ == "__main__":
	main()
