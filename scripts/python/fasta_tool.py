#!/usr/bin/env python
#==============================================================================
# NIBR
#
# @author: Michael Jones
# @date:   12/18/17 
#==============================================================================#
'''

 Describe

'''
__author__ = 'jonesmic'

from itertools import groupby
import re
import argparse

def fasta_iter(fasta_name):
    "first open the file outside "

    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)



def remove_reverse(fasta_file, out_fasta):
    it = fasta_iter(fasta_file)

    out_file  = open(out_fasta, 'w')

    for record in it:
        # Format FASTA for MASCOT - I think the TPP format may work for these.
        header = record[0]
        sequence = record[1]
        if 'Reverse' not in header:
            record_str = '>{0}\n{1}\n'.format(header, sequence)
            out_file.writelines(record_str)

    out_file.close()


def main(fasta_file, out_fasta):
    msg = 'Processing smORF FASTA files {0} to {1}'.format(fasta_file, out_fasta)
    print(msg)

    remove_reverse(fasta_file, out_fasta)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Remove the reverse sequences')
    parser.add_argument('--fasta_file', metavar='path', required=True,
                        help='The path to the FASTA file to parse. Specify multiple times if needed')
    parser.add_argument('--out_fasta', metavar='path', required=True,
                        help='path the created FASTA')

    args = parser.parse_args()
    main(fasta_file=args.fasta_file, out_fasta=args.out_fasta)

