
import urllib, xmltodict
import sys

'''
class for a uniprot feature containing:
feature type, description, location (start,end), uniprot id, uniprot name
'''
class Feature:
	def __init__(self, uniprot_id, t, description, location):
		self.uniprot_id = uniprot_id
		self.feature_type = t
		self.description = description
		self.location = [a-1 for a in location] #convert to 0-indexing
		self.start = int(self.location[0])
		self.end = int(self.location[1])
		self.length = self.end-self.start+1

class Uniprot_Entry:
	def __init__(self, uniprot_id, name, seq):
		self.id = uniprot_id
		self.name = name
		self.seq = seq
		self.features = []

	# for "in" comparisons
	def __eq__(self, ue_2):
		return self.id == ue2.id

	def add_feature(self, f):
		self.features.append(f)

	# helper methods
	def peptide_index(self, peptide_seq):
		peptide = ''.join([a for a in peptide_seq if a.isalpha()])
		return self.seq.find(peptide)

	def mod_index_in_peptide(self, peptide_seq_with_mod):
		return peptide_seq_with_mod.find("*")-1

	def mod_index_in_protein(self, 	peptide):
		return self.peptide_index(peptide)+self.mod_index_in_peptide(peptide)
	




'''
fetch sequenece and feature data from uniprot XML
returns string that is the protein seq, and a list of features as dicts
'''
def get_feature_data_from_uniprot_xml(uniprot_id):
	try :
		xml = urllib.urlopen("http://www.uniprot.org/uniprot/%s.xml" % uniprot_id.strip()).read()
		uniprot_data = xmltodict.parse(xml)['uniprot']['entry']
	except:
		return None
	
	# first get the sequence
	seq = ''.join([a.strip() for a in uniprot_data['sequence']['#text']])
	# then get name
	if 'recommendedName' in uniprot_data['protein']:
		name = uniprot_data['protein']['recommendedName']['fullName']
	else:
		if type(uniprot_data['protein']['submittedName']) == list:
			name = uniprot_data['protein']['submittedName'][0]['fullName']
		else:
			name = uniprot_data['protein']['submittedName']['fullName']
	if type(name) == xmltodict.OrderedDict:
		name = name['#text']
	# Create new uniprot object
	Uniprot = Uniprot_Entry(uniprot_id, name, seq)
	# now get all features
	if 'feature' in uniprot_data:
		u_features = uniprot_data['feature']
	else:
		return None

	if type(u_features) is not list:
		u_features = [u_features]

	# now iterate through each feature and get the relevant data
	for feature in u_features:
		f_type = feature['@type']
		try:
			f_description = feature['@description']
		except KeyError:
			f_description = ''
		location = feature['location']
		if 'begin' in location.keys():
			try:
				begin = int(location['begin']['@position'])
			except KeyError:
				begin = 0
			try:
				end = int(location['end']['@position'])
			except KeyError:
				end = len(seq) + 1         
			f_location = [begin,end]
		else:
			f_location = [int(location['position']['@position']),int(location['position']['@position'])]
		Uniprot.add_feature(
	    	Feature(uniprot_id, f_type, f_description, f_location)
	    	)

	return Uniprot

if __name__ == '__main__':
	analyzed_file = open(sys.argv[1])
	out_file = open(sys.argv[2], 'w')
	# throw out header
	analyzed_file.readline()
	uniprot_entries = {}
	feature_dict = {}
	stats = {}
	for line in analyzed_file:
		data = line.split(",")
		peptide = data[0].strip()
		area_ratio = float(data[2])
		uniprot_ids = set(data[3].split(' '))
		
		#out_file.write("\n%s\t%s\n" % (peptide, area_ratio))
		feature_dict[peptide] = (area_ratio,{})
		features = []
		for u_id in uniprot_ids:
			if u_id in uniprot_entries:
				uniprot = uniprot_entries[u_id]
				features = uniprot.features
			else:	
				uniprot = get_feature_data_from_uniprot_xml(u_id.strip())
				# if the uniprot xml can be fetched then add it to the dictionary
				if uniprot:
					uniprot_entries[u_id] = uniprot
					features = uniprot.features

			# if both classes were properly initialized
			if uniprot and features:
				for feature in features:
					mod_index = uniprot.mod_index_in_protein(peptide)
					if feature.start <= mod_index <= feature.end:
						if feature.length < 30:
							stats[feature.feature_type] = stats.setdefault(feature.feature_type, 0) + 1
						feature_dict[peptide][1].setdefault(u_id, []).append(feature)
						#out_file.write( "\t\t%s\t%s\t%s%s\t%s\t%s\n" % (feature.feature_type,feature.description,mod_index,feature.protein_seq[mod_index],feature.location,feature.end-feature.start+1))
	for ftype, count in stats.items():
		out_file.write("%s\t%s\n" % (ftype, count))
	out_file.write("\nPeptide\tArea Ratio\tUniprot ID\tAnnotation type\tAnnotation Description\tLocation\tAnnotation Span\tAnnotation Length\n")

	annotation_count = 0
	for peptide in feature_dict.keys():
		area_ratio = feature_dict[peptide][0] 
		out_file.write("%s\t%s\n" % (peptide, area_ratio))
		annotated = False 
		for u_id, features in feature_dict[peptide][1].items():
			uniprot = uniprot_entries[u_id]
			mod_index = uniprot.mod_index_in_protein(peptide)

			relevant_features = [f for f in features if f.length < 30]
			if len(relevant_features) > 0:
				annotated = True
				out_file.write("\t\t%s\t%s\n" % (u_id, uniprot_entries[u_id].name))
				for feature in relevant_features:
					out_file.write( "\t\t\t%s\t%s\t%s%s\t%s\t%s\n" % (feature.feature_type,feature.description,uniprot.seq[mod_index],mod_index+1,feature.location,feature.length))
			else:
				out_file.write("\t\t%s\t%s\tNo relevant annotation\t%s%s\n" % (u_id, uniprot_entries[u_id].name,uniprot.seq[mod_index],mod_index+1))
				
		if annotated:
			annotation_count += 1
		#out_file.write("\t\t\t\t\t\t\t\t%s\t%s\n" % (area_ratio, int(annotated)))


	out_file.write("\nTotal Peptides: %s\nAnnotated: %s\n\n" % (len(feature_dict.keys()), annotation_count))




			





