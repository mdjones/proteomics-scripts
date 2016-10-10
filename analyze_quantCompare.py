'''
Follows the standard analysis pathways from the QuantCompare program output file

1) read in peptides
	a) remove irrelevant modifications (Met Oxidation)
	b) remove probe modification, but store index in peptide
	c) keep track of how many runs it appeared in ("Unique column") and Area_ratio
2)	remove peptides in < 2 runs
3)	group peptides that contain one another, ex: ABC*DEFG and C*DEFG would be the same
4) 	average area rations and report: peptide, location of mod, area ratio, annotation (protein)
'''

import sys, re, argparse
from string import ascii_uppercase # simple list of uppercase letters

excepted_modifications = ["15.994915"]

'''
class for storing peptide data
contains close_match method to equate peptides that are truncated or extended
'''
class Peptide:
	def __init__(self, seq, mod_locs, ptm_is, area_ratio, annotation, run_counter):
		self.seq = seq
		self.mod_locs = [i-1 for i in mod_locs]
		self.area_ratio = area_ratio
		self.annotation = annotation
		self.uniprot_ids = []
		self.ptm_indices = ptm_is
		self.run_counter = run_counter
		# find all the valid uniprot ids that it matched with 
		for ann in annotation.split("#"):
			ann = ann.strip()
			match = re.match('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', ann) # regex for matching uniprot IDs from http://www.uniprot.org/help/accession_numbers
			if match:
				self.uniprot_ids.append(match.group())


	'''
	If one of the peptides is contained within another one, they are equivalent IF:
	the have diffmods in the same place
	'''
	def close_match(self, pep2):
		if pep2 == None:
			return False

		if pep2.seq in self.seq: # is pep2 within this peptide?
			i = self.seq.find(pep2.seq)
			if [a+i for a in pep2.mod_locs] == self.mod_locs: # are the diffmods the same (when aligned)
				return True 
		if self.seq in pep2.seq: # and vice versa
			i = pep2.seq.find(self.seq)
			if [a+i for a in self.mod_locs] == pep2.mod_locs:
				return True
		return False
		
	def offset(protein_seq):
		return protein_seq.find(self.seq)
	
	# is one or more of the mods spanned by a given feature which is given as a dict from the method
	# takes a peptide sequence to match to the peptide


	def __eq__(self, pep2):
		return self.close_match(pep2)
	def __ne__(self, pep2):
		return not self.close_match(pep2)


	def peptide_string(self):
		s = ""
		for i,c in enumerate(self.seq):
			s += c
			if i in self.mod_locs:
				s += "*"
		return s

	def uniprot_ids_str(self):
		return ' '.join(self.uniprot_ids)

	def run_count(self):
		return sum([int(a) for a in self.run_counter])




'''
takes in peptide, removes modifications, keeps track of where any mod
not in "excepted_modifications" is located
'''
def process_modifications(pep, ex_mods, ptm_ind = None):
	offset = None

	for ex in ex_mods:
		seq = pep.replace("(%s)"  % ex, "")


	#pep_split = re.split('[\(\)]', pep) # split on parenthsis
	# assume only one modified residue per peptide
	mod_locs = []
	i = 0
	for char in seq:
		if char == '(':
			mod_locs.append(i)
		if char in ascii_uppercase:
			i += 1
	seq = [p for p in seq if not p.isdigit() and p.isalpha()]
	peptide = "".join(seq)
	return peptide, mod_locs, offset

'''
fetch the data from the line and create a new peptide instance with relevant data
filter out peptides that appeared in only one run at the end, after merging
'''
def get_peptide_data_if_eligible(data):
	
	run_data = data[2].split(";")[:-1] # remove empty last spot
	
	run_counter = [a != 'X' for a in run_data]
	
	no_show_count = run_counter.count(False)
	
	total_runs = len(run_counter)
	#rc = total_runs - no_show_count
	#if no_show_count >= total_runs-1: # if it showed up in > 1 run - use it
	#	# 
	#	return None 
	
	ptm_idxs = [a for a in set(data[1].split("#")) if not a.startswith("M") and not a == "NA"]
	seq, mod_locs, offset = process_modifications(data[0], excepted_modifications)
	area_ratio  = float(data[5])
	annotation = data[22]
	peptide = Peptide(seq, mod_locs, ptm_idxs, area_ratio, annotation, run_counter)

	return peptide

'''
Adds the peptide to the big nested list of modified peptides
keeps track of all "identical" peptides for review

peptide is considered identical if it matches a current peptide with +/- 1 aa on an end - still TODO
'''
def add_peptide_to_dict(p_dict, pep):
	for pep_group in p_dict:
		if pep in pep_group: # "in" uses built-in __eq__ and __ne__ methods of pep to assess equivalence
			pep_group.append(pep)
			return p_dict
	p_dict.append([pep])
	return p_dict

def pep_group_run_count(pep_group):
	total_rcr = [False] * len(pep_group[0].run_counter)
	for pep in pep_group:
		rc = pep.run_counter
		for i, v in enumerate(rc):
			if v:
				total_rcr[i] = True
	total_rc = sum([int(a) for a in total_rcr])
	return total_rc



def out_line(peptide, ar):
	#if peptide.global_mod_locs:
	#	global_mod_locs = " ".join([str(a) for a in peptide.global_mod_locs])
	#	return "%s,%s,%s,%s,%s\n" % (peptide.peptide_string(), ar, global_mod_locs, peptide.uniprot_ids_str(), peptide.annotation)
	#else:
	return "%s,%s,%s,%s,%s,%s,%s\n" % (peptide.peptide_string(),' '.join(peptide.ptm_indices), ar, peptide.uniprot_ids_str(), peptide.annotation, peptide.run_counter, peptide.run_count())
'''
write a file with the desire data from each peptide
'''
def write_data_from_peptide_list(pep_list,outfile,reverse = False, verbose = False):
	outfile.write("%s,%s,%s,%s,%s\n" % ("Peptide", "PTM index from ip2", "average area ratio", "uniprot", "protein"))
	for pep_group in pep_list:
		run_count = pep_group_run_count(pep_group)

		pep = pep_group[0]
		peptide_seq = pep.seq
		mod_locs = pep.mod_locs
		annotation = pep.annotation
		ar = 0
		total_original_runs = 0
		for pep_a in pep_group:
			total_original_runs += pep_a.run_count()
			ar += pep_a.area_ratio * pep_a.run_count()
		ar = ar/(total_original_runs) # calculate average area ratio
		if reverse:
			try:
				ar = 1./ar
			except ZeroDivisionError:
				ar = 'inf'
		
		if verbose: # print all peptides from a group
			for p in pep_group:
				outfile.write("%s,%s,%s,%s,%s,%s\n" % (p.peptide_string(),' '.join(p.ptm_indices), str(ar) + "," + str(p.area_ratio), p.uniprot_ids_str(), p.annotation, p.run_count()))
		else:	
			outfile.write("%s,%s,%s,%s,%s,%s\n" % (pep.peptide_string(),' '.join(pep.ptm_indices), ar, pep.uniprot_ids_str(), pep.annotation, run_count))

'''
modifies the lines read in from the csv so as to be easily parseable
'''
def fix_csv_line(line):
	within_quotes = re.findall('\".*?\"',line) # find and capture all the blocks between quotes
	for block in within_quotes:
		line = line.replace(block, block.replace(",","#")) # and replace the commas with #s so comma delimiting works
	return line

def parse_input():
	parser = argparse.ArgumentParser()
	parser.add_argument("in_file", help = 'Quant or QuantCompare file to read')
	parser.add_argument("out_file", help = 'Name of csv file to write')
	parser.add_argument("-v","--verbose", help = 'print all peptides considered to be equivalent', action='store_true')
	parser.add_argument("-r","--reverse", help = 'reciprocal of area ratio', action='store_true')

	args = parser.parse_args()

	return open(args.in_file), open(args.out_file, 'w'), args.reverse, args.verbose



if __name__ =='__main__':
	quantc_file, out_file, reverse, verbose = parse_input()
	#quantc_file.readline() #throw out header line

	pep_list = []

	for line in quantc_file:
		'''
		Lots of clean up of the input text, replaces internal commas in fields with #
		'''
		line = fix_csv_line(line).strip()
		data = line.split(",") # split on " because of terrible formatting by ip2 and no XML availability
		data = [a.strip("\"") for a in data] # strip off excess quotes
		

		peptide = get_peptide_data_if_eligible(data) # create a new peptide instance if the peptide appears in > 1 run
		if peptide != None:
			add_peptide_to_dict(pep_list, peptide) # add it to the list of peptides
		

	# cut off all peptides that appeared in < 2 runs
	pep_list = [a for a in pep_list if pep_group_run_count(a) >= 2]
	write_data_from_peptide_list(pep_list, out_file, reverse, verbose) # write the list of peptides to file




