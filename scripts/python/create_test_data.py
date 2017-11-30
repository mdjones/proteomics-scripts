#!/usr/bin/env python

__author__ = 'JONESMIC'

#!/usr/bin/env python
# ==============================================================================
#
#
# @author: Michael Jones
# @date:   11/30/17�
# ==============================================================================#
'''

 Create test data from real data. 

'''

import re
import random
import argparse
import pandas as pd
import os
import subprocess

ip2_pept_pattern = re.compile(r'(\w+)(C\(\d+\.\d+\))(\w+)')

def scramble_string(word):
    word = list(word)
    random.shuffle(word)
    return ''.join(word)

def scramble_ip2(ip2_pep):
    match = re.match(ip2_pept_pattern, ip2_pep)
    scrambled_pep = '{0}{1}{2}'.format(scramble_string(match.group(1)), match.group(2), scramble_string(match.group(3)))
    return scrambled_pep

def scramble_protein(protein):
    proteins = protein.split(',')

    result = []

    for protein in proteins:
        scrambled_prot = []
        if 'Reverse_tr' in protein:
            protein = protein.replace('Reverse_tr|', '')
            scrambled_prot.append('Reverse_tr|')

        protein = protein.split()
        for word in protein:
            #scrambled_prot.append(scramble_string(word))
            scrambled_prot.append(word)
        result.append(' '.join(scrambled_prot))

    return ', '.join(result)

def scramble_data(peptide_list_file, file_type):
    df = None

    sequence_col = None
    protein_col = None
    if file_type == 'CSV':
        df = pd.read_csv(peptide_list_file)
        sequence_col = 'sequence'
        protein_col = 'protein'
    elif file_type == 'TXT':
        df = pd.read_table(peptide_list_file, skiprows=24)
        sequence_col = 'SEQUENCE'
        protein_col = 'PROTEIN'




    #df[sequence_col] = df[sequence_col].apply(scramble_ip2)
    #df[protein_col] = df[protein_col].apply(scramble_protein)



    return df

def extract_by_sequence_pattern(out_dir, peptide_list_file, sequence_patterns_file):
    #grep -E -f ${DIR}/peptideListGrepPatterns.txt ${PEP_LIST_IN} > ${TEST_DATA_DIR}/peptideList.csv

    out_file = '{0}/{1}'.format(out_dir, "peptideList.csv")

    f = open(out_file, 'w')
    cmd = ['egrep', '-f', sequence_patterns_file, peptide_list_file]
    print(subprocess.run(cmd, stdout=f))


def process_with_analyzeQuantComapare(out_dir):

    python_bin = '/usr/bin/python'
    analyze_quant_comapare = './scripts/nomurarg-scripts/analyze_quantCompare.py'

    in_file = '{0}/{1}'.format(out_dir, "peptideList.csv")
    out_file = '{0}/{1}'.format(out_dir, "results.csv")
    ## "${PYTHON} ${ANALYZE_QUANT_SCRIPT} ${TEST_DATA_DIR}/peptideList.csv ${TEST_DATA_DIR}/results.csv"
    cmd = [python_bin, analyze_quant_comapare, in_file, out_file]
    print(subprocess.check_output(cmd))

    out_file = '{0}/{1}'.format(out_dir, "resultsverbose.csv")
    ## "${PYTHON} ${ANALYZE_QUANT_SCRIPT} -v ${TEST_DATA_DIR}/peptideList.csv ${TEST_DATA_DIR}/resultsverbose.csv"
    cmd = [python_bin, analyze_quant_comapare, '-v', in_file, out_file]
    print(subprocess.check_output(cmd))


def main(out_dir, peptide_list_file, sequence_patterns_file):
    file_type = None
    if peptide_list_file.endswith('csv'):
        file_type = 'CSV'
    elif peptide_list_file.endswith('csv'):
        file_type = 'TXT'

    extract_by_sequence_pattern(out_dir, peptide_list_file, sequence_patterns_file)

    testPeptideListFile = '{0}/{1}'.format(out_dir, "peptideList.csv")
    df = scramble_data(testPeptideListFile, file_type)

    ## Write file
    if file_type == 'CSV':
        df.to_csv(testPeptideListFile, index=False)
    elif file_type == 'TXT':
        df.to_csv(testPeptideListFile, index=False, sep="\t")

    process_with_analyzeQuantComapare(out_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create test data.')

    parser.add_argument('--out_dir', type=str, required=True,
                        help='The directory to put the test data into')

    parser.add_argument('--peptide_list_file', type=str, required=True,
                        help='A peptideList.csv AKA quant_stat_peptide_compare.txt')

    parser.add_argument('--sequence_patterns', type=str, required=True,
                        help='A file listing the IP2 format pep patterns. Example: "ICD\([0-9]+.[0-9]+\)LLLLLEK"')

    args = parser.parse_args()

    main(args.out_dir, args.peptide_list_file, args.sequence_patterns)

