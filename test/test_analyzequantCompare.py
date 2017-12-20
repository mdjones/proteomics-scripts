#!/usr/bin/env python
#==============================================================================
# NIBR
#
# @author: Michael Jones
# @date:   5/23/17 
#==============================================================================#
'''

 Unit tests to test the meta-data package creation

'''
__author__ = 'jonesmic'

import  numpy as np
from nbcpact.ucbre import AnalyzeQuantCompare, UcbreUtils, DataAccessObject, PeptidesFromPeptideListBuilder, Peptide, PeptideGroup

import pkg_resources
from nose.tools import nottest
import warnings

USE_BIG_TEST = False


def setup():
    print ("SETUP! Validate")

def teardown():
    print ("TEAR DOWN! Validate")

def test_basic():
    print ("I ran a test in Validate")


def get_file_paths():

    path_dict = {}
    resource_package = __name__  # Could be any module/package name

    suffix = '_EN80.csv' if USE_BIG_TEST else '.csv'
    resource_path = 'data/peptideList{0}'.format(suffix)
    path_dict['peptideList'] = pkg_resources.resource_filename(resource_package, resource_path)

    resource_path = 'data/results{0}'.format(suffix)
    path_dict['results'] = pkg_resources.resource_filename(resource_package, resource_path)

    resource_path = 'data/resultsverbose{0}'.format(suffix)
    path_dict['resultsverbose'] = pkg_resources.resource_filename(resource_package, resource_path)

    return path_dict

def test_results_data_creation():
    path_dict = get_file_paths()

    peptide_generator = PeptidesFromPeptideListBuilder(peptide_list_file=path_dict['peptideList'])
    analyzeQuantCompare = AnalyzeQuantCompare(peptide_generator=peptide_generator)

    groups = analyzeQuantCompare.build_peptide_groups()
    generated_results_csvDF = AnalyzeQuantCompare.build_results_from_peptide_groups(groups)

    ucbResultsDF = UcbreUtils.read_results_csv(path_dict['results'])

    compare_dataframes(ucbdf=ucbResultsDF, novdf=generated_results_csvDF,
                       merge_cols=['Peptide'],
                       identical_cols=['run_count', 'uniprot', 'annotations'],
                       close_cols=['mean_group_ratio'],
                       warn_cols=['ptm_index_from_ip2'])


def test_results_data_creation_from_peptide_generator():
    path_dict = get_file_paths()
    peptide_generator = PeptidesFromPeptideListBuilder(peptide_list_file=path_dict['peptideList'])
    analyzeQuantCompare = AnalyzeQuantCompare(peptide_generator=peptide_generator)

    groups = analyzeQuantCompare.build_peptide_groups()
    generated_results_csvDF = AnalyzeQuantCompare.build_results_from_peptide_groups(groups)

    ucbResultsDF = UcbreUtils.read_results_csv(path_dict['results'])

    compare_dataframes(ucbdf=ucbResultsDF, novdf=generated_results_csvDF,
                       merge_cols=['Peptide'],
                       identical_cols=['run_count', 'uniprot', 'annotations'],
                       close_cols=['mean_group_ratio'],
                       warn_cols=['ptm_index_from_ip2'])

def test_results_verbose_data_creation():
    path_dict = get_file_paths()

    peptide_generator = PeptidesFromPeptideListBuilder(peptide_list_file=path_dict['peptideList'])
    analyzeQuantCompare = AnalyzeQuantCompare(peptide_generator=peptide_generator)

    groups = analyzeQuantCompare.build_peptide_groups()
    generated_results_csvDF = AnalyzeQuantCompare.build_results_from_peptide_groups(groups, verbose=True)

    ucbResultsDF = UcbreUtils.read_results_verbose_csv(path_dict['resultsverbose'])

    compare_dataframes(ucbdf=ucbResultsDF, novdf=generated_results_csvDF,
                       merge_cols=['Peptide', 'ratios'],
                       identical_cols=['run_count', 'uniprot', 'annotations'],
                       close_cols=['mean_group_ratio'],
                       warn_cols=['ptm_index_from_ip2'])

## Helpers
def compare_dataframes(ucbdf=None, novdf=None, merge_cols=None, identical_cols=None, close_cols=None,
                           warn_cols=None):
        """
        Throw an error if any values are not correct.
        :param ucbdf:
        :param novdf:
        :param merge_cols:
        :param identical_cols: columns that must have the same exact value.
        :param close_cols: columns that are the same within the default np.close tolerance.
        :param warn_cols: columns that should be close but there might be a bug in the UBC code.
        """

        assert ucbdf.index.size == novdf.index.size, \
            'Data frames do not have the same number of rows ucbdf={0} and novdf={1}'.format(ucbdf.index.size,
                                                                                             novdf.index.size)

        mergedDF = ucbdf.merge(novdf, on=merge_cols, suffixes=('_ucb', '_nov'), how='outer')

        for col in identical_cols:
            nov_col = '{0}_nov'.format(col)
            ubc_col = '{0}_ucb'.format(col)
            diff_df = mergedDF[mergedDF[nov_col] != mergedDF[ubc_col]]
            assert diff_df.empty, 'Non matching values for {0} --- {1}'.format(col, diff_df.head())

        for col in close_cols:
            nov_col = '{0}_nov'.format(col)
            ubc_col = '{0}_ucb'.format(col)
            diff_df = mergedDF[~(np.isclose(mergedDF[nov_col], mergedDF[ubc_col]))]
            assert diff_df.empty, 'Non matching values for {0} --- {1}'.format(col, diff_df.head())

        for col in warn_cols:
            nov_col = '{0}_nov'.format(col)
            ubc_col = '{0}_ucb'.format(col)
            diff_df = mergedDF[mergedDF[nov_col] != mergedDF[ubc_col]]
            if not diff_df.empty:
                warnings.warn('Non matching values for {0} --- {1}'.format(col, diff_df.head()))

def test_peptide_equality():
    """
     The peptide.close_match function may not give expected results in each case.
     Example: It may fail on this partial overlap.

             ABPDEKGHIJcK
                   GHIJcK
                   GHIJcKLMNOPQR

    """

    peptides = ['ABPDEKGHIJcK', 'GHIJcK', 'GHIJcKLMNOPQR']

    print('Input peptides {0}'.format(peptides))
    peptide_generator = MockPeptideGenerator(peptide_mod_seqs = peptides)
    do_grouping(peptide_generator)

    print('************************************************************')
    peptides = ['ABPDEKGHIJcK', 'GHIJcKLMNOPQR']
    print('Input peptides {0}'.format(peptides))
    peptide_generator = MockPeptideGenerator(peptide_mod_seqs=peptides)
    do_grouping(peptide_generator)



def do_grouping(peptide_generator):
    peptides = peptide_generator.generate_peptides()

    peptide_groups = []

    sequences = []
    for pep in peptides:
        sequences.append(pep.sequence)
        for pep_group in peptide_groups:
            if pep_group.contains_peptide(
                    pep):  # "in" uses built-in __eq__ and __ne__ methods of pep to assess equivalence
                pep_group.append(pep)

        ## Not found so create new peptide group
        peptideGroup = PeptideGroup(pep)
        peptide_groups.append(peptideGroup)


    group_num = 0
    for group in peptide_groups:
        group_peps = []
        group_num += 1
        for peptide in group.peptides():
            group_peps.append(peptide.sequence)
        print('Group{} {}'.format(group_num, group_peps))

class MockPeptideGenerator:

    def __init__(self, peptide_mod_seqs=None):
        """

        :param peptide_mod_seqs: a list of sequences with the Cys modifcation site in lower case. ['ABPDEKGHIJcK', 'GHIJcK', 'GHIJcKLMNOPQR']
        """
        self.__peptide_mod_seqs = peptide_mod_seqs

    def generate_peptides(self):

        peptides = []

        for mod_seq in self.__peptide_mod_seqs:
            peptide = Peptide(sequence=mod_seq.upper(),
                              mod_locs=[mod_seq.index('c')],
                              ptm_indices=None,
                              area_ratio=1.3,
                              area_ratios=[1.2, 1.4, 1.3],
                              annotation=None,
                              uniprot_ids=None,
                              run_counter=None,
                              decoy=None,
                              unique1=None,
                              ip2_peptide=None)

            peptides.append(peptide)

        return peptides


