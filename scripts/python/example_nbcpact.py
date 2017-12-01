#!/usr/bin/env python
#==============================================================================
# NIBR
#
# @author: Michael Jones
# @date:   12/1/17 
#==============================================================================#
'''

 A script to encapsulate basic usage of the library

'''
__author__ = 'jonesmic'

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='nbcpact-scripts.')

    parser.add_argument('--type', type=str, required=True,
                        help='The directory to put the test data into')

    parser.add_argument('--out_dir', type=str, required=True,
                        help='The directory to put the test data into')

    parser.add_argument('--peptide_list_file', type=str, required=True,
                        help='A peptideList.csv AKA quant_stat_peptide_compare.txt')

    parser.add_argument('--sequence_patterns', type=str, required=True,
                        help='A file listing the IP2 format pep patterns. Example: "ICD\([0-9]+.[0-9]+\)LLLLLEK"')

    args = parser.parse_args()

    main(args.out_dir, args.peptide_list_file, args.sequence_patterns)

