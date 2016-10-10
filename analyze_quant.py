'''
Same as analyze_quantCompare but uses just the quantification from a single run

Takes the PTM index from the file and carries it through

only how the peptide is fetched from the data is changed

peptide class natively has support for PTM locations for further implementation
'''
import sys, re

# import all relevant classes and methods from other script
# makes changing things easier	
from analyze_quantCompare import Peptide
from analyze_quantCompare import excepted_modifications
from analyze_quantCompare import process_modifications, add_peptide_to_dict, out_line, write_data_from_peptide_list, fix_csv_line, parse_input



'''
fetch the data from the line and create a new peptide instance with relevant data integrate the PTM location
into the peptide object
'''
def get_peptide_data_from_quant(data):
	# split every ptm index out using re, IP2 collects non unquie and uses semicolons and commas as separators
	# re is SO SLOW 
	
	seq, mod_locs, offset = process_modifications(data[1].strip("singleton"), excepted_modifications, re.split("\W+",data[2].strip()))
	area_ratio  = float(data[14])
	annotation = data[17].strip("[").strip("]")
	peptide = Peptide(seq, mod_locs, area_ratio, annotation)
	return peptide



if __name__ =='__main__':
	quant_file, out_file, reverse, verbose = parse_input()
	quant_file.readline() #throw out header line

	pep_list = []

	for line in quant_file:
		'''
		Lots of clean up of the input text
		'''
		line = fix_csv_line(line).strip()
		data = line.split(",") # split on " because of terrible formatting by ip2 and no XML availability
		data = [a.strip("\"") for a in data] # strip off excess quotes
		

		peptide = get_peptide_data_from_quant(data) # create a new peptide instance if the peptide appears in > 1 run

		if peptide != None:
			add_peptide_to_dict(pep_list, peptide) # add it to the list of peptides
		


	write_data_from_peptide_list(pep_list, out_file, reverse, verbose) # write the list of peptides to file







