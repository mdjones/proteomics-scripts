#!/usr/bin/env python
# ==============================================================================
# NIBR
#
# @author: Michael Jones
# @date:   11/14/17 
# ==============================================================================#
'''

 A set of tools to reverse engineer or replicate ucb analysis

'''
__author__ = 'jonesmic'
import re
from string import ascii_uppercase

import numpy as np
import pandas as pd


class Peptide:
    """
        This is adopted from from analyze_quantCompare.py
    """
    annotation = None
    area_ratios = None
    sequence = None
    uniprot_ids = None
    ptm_indices = None
    mod_locs = None
    ip2_peptide = None
    unique1=None
    run_counter=None

    def __init__(self, ip2_peptide=None, sequence=None, mod_locs=None, ptm_indices=None, area_ratio=None,
                 area_ratios=None, annotation=None, uniprot_ids=None, run_counter=None, decoy=None,
                 unique1=None):
        """
        Create a peptide from the peptideList file. Each object represents one row in the file.
        :param ip2_peptide: Ex:  AFAFVTFADDQIAQSLC(470.29977)GEDLIIK this is the primary grouping key
        :param sequence: Ex:  AFAFVTFADDQIAQSLCGEDLIIK
        :param mod_locs: Ex:  [17]
        :param ptm_indices: Ex:  ['C244;C244;C244;C244;C244']
        :param area_ratio: Ex:  1.78394 -- This comes from the AREA_MEDIAN_RATIO_1 in the peptideList.csv. It also is the median of the ratios.
        :param area_ratios: Ex:  [3.2027, 1.78394, 1.76109]
        :param annotation: Ex:  A0A087X260 TAR DNA-binding protein 43 OS=Homo sapiens GN=TARDBP PE=1 SV=1 , A0A087WYY0 TAR DNA-binding protein 43 OS=Homo sapiens GN=TARDBP PE=1 SV=1 , B1AKP7 TAR DNA-binding protein 43 OS=Homo sapiens GN=TARDBP PE=1 SV=1 , Q13148 TAR DNA-binding protein 43 OS=Homo sapiens GN=TARDBP PE=1 SV=1 , G3V162 TAR DNA binding protein, isoform CRA_d OS=Homo sapiens GN=TARDBP PE=1 SV=1 ,
        :param uniprot_ids: Ex:  ['A0A087X260', 'A0A087WYY0', 'B1AKP7', 'Q13148', 'G3V162']
        :param run_counter: Ex:  [True, True]
        :param decoy: Ex:  False
        :param unique1: Ex:  ,,;,;
        """
        self.ip2_peptide = ip2_peptide
        self.sequence = sequence
        # self.mod_locs = [i - 1 for i in mod_locs]
        self.mod_locs = mod_locs
        self.area_ratios = area_ratios
        self.annotation = annotation
        self.uniprot_ids = uniprot_ids
        self.ptm_indices = ptm_indices
        self.run_counter = run_counter
        self.decoy = decoy
        self.area_ratio = area_ratio
        self.unique1 = unique1

        self.__validate()

    '''
    If one of the peptides is contained within another one, they are equivalent IF:
    the have diffmods in the same place
    '''

    def close_match(self, pep2):
        if pep2 == None:
            return False

        if pep2.sequence in self.sequence:  # is pep2 within this peptide?
            i = self.sequence.find(pep2.sequence)
            if [a + i for a in pep2.mod_locs] == self.mod_locs:  # are the diffmods the same (when aligned)
                return True
        if self.sequence in pep2.sequence:  # and vice versa
            i = pep2.sequence.find(self.sequence)
            if [a + i for a in self.mod_locs] == pep2.mod_locs:
                return True
        return False

    def __eq__(self, pep2):
        return self.close_match(pep2)

    def __ne__(self, pep2):
        return not self.close_match(pep2)

    def peptide_string(self):
        peptide = self.sequence
        i = 0
        for mod_loc in self.mod_locs:
            peptide = peptide[:mod_loc+i] + '*' + peptide[mod_loc+i:]
            i += 1
        return peptide


    def uniprot_ids_str(self):
        return ' '.join(self.uniprot_ids)

    def run_count(self):
        return sum(self.run_counter)

    def __validate(self):
        calc_median = np.median(np.array(self.area_ratios))
        assert np.isclose(self.area_ratio, calc_median), \
            'Area ratios do not represent the reported median ratio {0}, area_ratio {1}, area_ratios {2}'.format(
                self.ip2_peptide, self.area_ratio, self.area_ratios)


class PeptideGroup:

    def __init__(self, founder_pep):
        self.__peptides = [founder_pep]
        #self.__group_run_count = self.__pep_group_run_count()

    def append(self, peptide):
        self.__peptides.append(peptide)
        #self.__group_run_count = self.__pep_group_run_count()

    def contains_peptide(self, peptide):
        return peptide in self.__peptides

    def peptides(self):
        return self.__peptides

    #def group_run_count(self):
    #    return self.__group_run_count

    def total_original_runs(self):
        total_original_runs = 0
        for pep in self.__peptides:
            total_original_runs += pep.run_count()
        return total_original_runs

    def all_ratios_in_group(self):
        ars = []
        for pep in self.__peptides:
            ars += pep.area_ratios
        return ars

    def pep_group_run_count(self):
        """
        Count the run if any peptide in this group was observed at least one time.


        Original algorithm from UCB

        def __pep_group_run_count(self):
            total_rcr = [False] * max([len(p.run_counter) for p in self.__peptides])
            for pep in self.__peptides:
                rc = pep.run_counter
                for i, v in enumerate(rc):
                    if v:
                        total_rcr[i] = True
            total_rc = sum([int(a) for a in total_rcr])
            return total_rc


        :return: number of runs any peptide in this group was observed.
        """

        total_rcr = [False] * max([len(p.run_counter) for p in self.__peptides])
        total_rc = np.nan
        for pep in self.__peptides:
            rc = pep.run_counter
            for i, v in enumerate(rc):
                if v:
                    total_rcr[i] = True
            total_rc = sum([int(a) for a in total_rcr])

        return total_rc

class UcbreUtils:

    @staticmethod
    def read_results_csv(file_path):
        """
        Read the CSV file and convert to a dataframe

        results csv is missing the last column header
        :param file_path:
        :return: a dataframe
        """

        header = ['Peptide', 'ptm_index_from_ip2', 'mean_group_ratio', 'uniprot', 'annotations', 'run_count']

        df = pd.read_csv(file_path, names=header, header=None, skiprows=1)
        df['annotations'] = df['annotations'].str.replace('#', ',')

        return df

    @staticmethod
    def read_results_verbose_csv(file_path):
        """
        Read the CSV file and convert to a dataframe

        results csv is missing the last column header
        :param file_path:
        :return: a dataframe
        """

        header = ['Peptide', 'ptm_index_from_ip2', 'mean_group_ratio', 'ratios', 'uniprot', 'annotations', 'run_count']

        df = pd.read_csv(file_path, names=header, header=None, skiprows=1)
        df['annotations'] = df['annotations'].str.replace('#', ',')

        return df

class AnalyzeQuantCompare:

    def __init__(self, peptide_list_file):
        """
        Initialize this analyis with a peptide list
        :param peptide_list_file:
        """
        self.__peptide_list_file = peptide_list_file

    @staticmethod
    def build_results_from_peptide_groups(groups, verbose=False, printextracols=False):
        """
        Build a DF representing the UCB results file format.
        :param verbose: if true print all peptides in the group
        :param printextra: If true print extra diagnosis columns
        :return: a datarame
        """
        columns = ['Peptide', 'ptm_index_from_ip2', 'mean_group_ratio', 'uniprot', 'annotations',
                   'run_count']

        verbose_columns = ['Peptide', 'ptm_index_from_ip2', 'mean_group_ratio', 'ratios',
                           'uniprot', 'annotations', 'run_count']


        data = {'ip2_sequence': [], 'sequence': [], 'Peptide': [], 'uniprot': [], 'annotations': [],
                'ptm_index_from_ip2': [],
                'run_count': [], 'all_ratios_in_group': [], 'ratios': [],
                'mean_group_ratio': [], 'mod_locs': [], 'unique1': []}


        for pep_group in groups:
            
            peps = pep_group.peptides() if verbose else [pep_group.peptides()[0]]
            
            for pep in peps:
                data['ip2_sequence'].append(pep.ip2_peptide)
                data['sequence'].append(pep.sequence)
                data['Peptide'].append(pep.peptide_string())
                data['annotations'].append(pep.annotation)
                data['uniprot'].append(pep.uniprot_ids_str())

                if verbose:
                    data['run_count'].append(pep.run_count())
                else:
                    data['run_count'].append(pep_group.pep_group_run_count())


                data['ptm_index_from_ip2'].append(pep.ptm_indices)
                data['mod_locs'].append(pep.mod_locs)
                data['unique1'].append(pep.unique1)
    
                all_ratios_in_group = pep_group.all_ratios_in_group()
                data['all_ratios_in_group'].append(all_ratios_in_group)

                data['mean_group_ratio'].append(np.mean(np.array(all_ratios_in_group)))
                data['ratios'].append(pep.area_ratios)

        df = pd.DataFrame(data=data)
        df['ptm_index_from_ip2'] = df['ptm_index_from_ip2'].apply(lambda x: ';'.join(x))

        if verbose:
            df['ratios'] = df['ratios'].apply(lambda x : ";".join(map(str, x)))

        if not printextracols:
            if verbose:
                df = df[verbose_columns]
            else:
                df = df[columns]
        return df



    def build_peptide_groups(self, minruncount=2):

        """
        From a peptideList.csv file extract all peptides.
        Filter decoy peptides
        Filter pep_group_run_count() >= minruncount
        :param minruncount: group must represent at least this number of min runs. Defaults to 2
        :return: a list of peptide groups
        """
        peptides = self.convert_peptideslist2peptides()

        filt_peps = list(filter(lambda p: not p.decoy, peptides))

        peptide_groups = []
        for peptide in filt_peps:
            self.__add_peptide_to_group(peptide_groups, peptide)

        peptide_groups = list(filter(lambda g: g.pep_group_run_count() >= minruncount, peptide_groups))

        return peptide_groups


    def __add_peptide_to_group(self, peptide_groups, pep):
        for pep_group in peptide_groups:
            if pep_group.contains_peptide(
                    pep):  # "in" uses built-in __eq__ and __ne__ methods of pep to assess equivalence
                pep_group.append(pep)
                return peptide_groups

        ## Not found so create new peptide group
        peptideGroup = PeptideGroup(pep)
        peptide_groups.append(peptideGroup)
        return peptide_groups


    def convert_peptideslist2peptides(self):
        """

        :param file_path: the location of the peptide list file
        :return: A list of Peptide objects for each line in the file.
        """

        df = self.read_peptide_list_file()
        peptides = self.__get_peptides_from_data_frame(df)
        return peptides

    def read_peptide_list_file(self):
        """
        Read the CSV file and convert to a dataframe

        Now that it has a header. Nothing really needs to be done.
        :param file_path:
        :return: a dataframe
        """
        df = pd.read_csv(self.__peptide_list_file)

        return df

    def __init_peptide(self, row):
        peptide = Peptide(sequence=row['peptide'],
                          mod_locs=row['mod_locs'],
                          ptm_indices=row['ptm_indices'],
                          area_ratio=row['file_AREA_MEDIAN_RATIO_1'],
                          area_ratios=row['area_ratios'],
                          annotation=row['annotation'],
                          uniprot_ids=row['uniprot_ids'],
                          run_counter=row['run_counter'],
                          decoy=row['decoy'],
                          unique1=row['UNIQUE_1'],
                          ip2_peptide=row['ip2_peptide'])

        return peptide

    '''
    takes in peptide, removes modifications, keeps track of where any mod
    not in "excepted_modifications" is located
    '''

    def __process_modifications(self, ip2_pep):
        ## TODO: maybe should be inclusive instead of exclusive. i.e. retain IsoTop Only
        excepted_modifications = ["15.994915"]

        for ex in excepted_modifications:
            ip2_pep = ip2_pep.replace("(%s)" % ex, "")

        # pep_split = re.split('[\(\)]', pep) # split on parenthsis
        # assume only one modified residue per peptide
        mod_locs = []
        i = 0
        for char in ip2_pep:
            if char == '(':
                mod_locs.append(i)
            if char in ascii_uppercase:
                i += 1

        return mod_locs

    def __get_peptides_from_data_frame(self, df):
        """
                Attempt to recreate Nomura / Karl analyze_quantCompare.py
                inputs a peptideList data frame see: dpro.data.get_ucb_peptide_list

                :param sequence: a pepideList Data Frame.
            """

        data = {}

        ## Clean up dataframe
        df = df[~df['PTM_INDEX'].isnull()]

        ## Transfer some values directly
        data['annotation'] = df['protein']
        ## TODO: Comapare this to the value calculated in Peptide
        data['file_AREA_MEDIAN_RATIO_1'] = df['AREA_MEDIAN_RATIO_1']
        ## TODO: Use to compare the run count
        data['UNIQUE_1'] = df['UNIQUE_1']

        # Clean up df
        ## Nothing to see here yet.

        # declare decoys
        data['decoy'] = (df['protein'].str.contains('Reverse_'))

        # Extract uniprot IDs from annotation
        """
        From http://www.uniprot.org/help/accession_numbers
        #'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})
        
        """
        uniprot_rx = re.compile(r'\b([OPQ]\d[A-Z0-9]{3}\d|[A-NR-Z]\d([A-Z][A-Z0-9]{2}\d){1,2})\b')
        data['uniprot_ids'] = list(df['protein'].str.extractall(uniprot_rx).iloc[:, 0].unstack().values)


        # Process mods
        data['mod_locs'] = df['sequence'].apply(self.__process_modifications)
        data['peptide'] = df['sequence'].str.replace('\([\d.]+\)', '')
        data['ip2_peptide'] = df['sequence']

        # Global mod indexes
        ## ptm_idxs = [a for a in set(data[1].split("#")) if not a.startswith("M") and not a == "NA"]
        ## Then later ' '.join(p.ptm_indices)
        data['ptm_indices'] = df['PTM_INDEX'].apply(lambda x: [a for a in x.split(',') if not a.startswith("M") and not a == "NA"])

        #data['ptm_indices'] = data['ptm_indices'].apply(lambda x : ' '.join(x))


        # TODO: Need to understand what the semicolon means.
        # Looks like area_ratio this is already in the peptide file. Why is it being recalulated?

        data['area_ratios'] = df['AREA_RATIO_ALL_1'].apply(lambda x: [float(a) for a in re.split('[,;]', x) if a and a != 'X'])

        ## Get area ratios withough the Xs
        data['run_data'] = df['UNIQUE_1'].apply(lambda x: x.split(";")[:-1])
        data['run_counter'] = data['run_data'].apply(lambda x: [a != 'X' for a in x])

        resultDF = pd.DataFrame(data)
        ## Remove None from Uniprot
        resultDF['uniprot_ids'] = resultDF['uniprot_ids'].apply(lambda L: [x for x in L if x is not None])


        peptides = resultDF.apply(self.__init_peptide, axis=1)


        return peptides


