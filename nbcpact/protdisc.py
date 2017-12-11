#!/usr/bin/env python
# ==============================================================================
# NIBR
#
# @author: Michael Jones
# @date:   12/8/17 
# ==============================================================================#
'''

 Library classes for accessing PD.THis should be used to a common lib once it is more mature.

'''
__author__ = 'jonesmic'

from enum import Enum
import re
import numpy as np

from sqlalchemy import create_engine, inspect
import pandas as pd
import struct


class DataType(Enum):
    Float = 1
    Integer = 2


class PDReader:
    __peptide_group_mod_pattern = re.compile(r'\[(\w\d+)\]')

    def __init__(self, pd_result_file=None,
                 num_quant_channels=2,
                 pd_version='2.1'):

        self.__num_quant_channels = num_quant_channels

        url = 'sqlite:///{0}'.format(pd_result_file)
        self.__engine = create_engine(url)

        self.__pd_version = pd_version
        self.__data_cache = {}

    def __extract_values(self, binary_data, n=None, dataType=None):
        """
        Adapted from Tomas Reijtar Unicorn project

        values: bytes from the blob
        n: number of channels
        t: type of data
            'decimal' for values in decimal format such as Abundances
            'integer_number' for values such as 'Found in'

        """

        result = []

        if binary_data:
            if dataType == DataType.Float:
                for i in range(n):
                    sub = binary_data[9 * i:9 * i + 8]
                    result.append(struct.unpack("d", sub)[0])
            else:
                for i in range(n):
                    result.append(struct.unpack("i", binary_data[5 * i:5 * i + 4]))

            if len(result) == 1:
                result = result[0]

        # take care of the missing values
        else:
            if dataType == DataType.Float:
                result = [0.0] * n
            elif (dataType == DataType.Integer):
                result = [0] * n

        return result

    def __get_found_raw_files(self, targetPeptideGroupsPeptideGroupID):
        df = self.get_target_psms()
        df = df[df['TargetPeptideGroupsPeptideGroupID'] == targetPeptideGroupsPeptideGroupID]

        values = df['SpectrumFileName'].values

        return values

    def __get_local_mod_locations(self, modifications):
        mod_locs = re.findall(self.__peptide_group_mod_pattern, modifications)
        return mod_locs

    def __get_global_mod_locations(self, row):
        local_mods = row['local_mod_locs']
        local_mods = [re.match(r'(\w)(\d+)', mod).groups() for mod in local_mods]
        peptideGroupID = row['PeptideGroupID']
        peptide_sequence = row['Sequence']

        df = self.__get_target_proteins()
        proteinSequences = df[df['TargetPeptideGroupsPeptideGroupID'] == peptideGroupID]['Sequence']

        peptide_starts = []
        for proteinSequence in proteinSequences.values:
            peptide_starts.append((proteinSequence.find(peptide_sequence)))

        global_positions = []
        for peptide_start in peptide_starts:
            pos_strings = []



            for local_mod in local_mods:
                local_position = int(local_mod[1])
                global_pos = np.nan if peptide_start == -1 else (local_position + peptide_start)
                pos_strings.append('{0}{1}'.format(local_mod[0],global_pos))

            global_positions.append(','.join(pos_strings))

        return global_positions

    def __read_data_frame(self, sql_string):
        connection = self.__engine.connect()
        df = pd.read_sql(sql_string, connection)
        connection.close()
        return df

    def __get_target_proteins(self):
        data_name = 'target_proteins'

        if not data_name in self.__data_cache.keys():
            sqlStr = """
                        SELECT
                        t1.TargetPeptideGroupsPeptideGroupID,Description,Sequence,Accession
                        FROM TargetPeptideGroupsTargetProteins t1,TargetProteins t2
                        WHERE t1.TargetProteinsUniqueSequenceID = t2.UniqueSequenceID
                    """
            self.__data_cache[data_name] = self.__read_data_frame(sqlStr)

        return self.__data_cache[data_name]

    def get_target_psms(self):
        data_name = 'target_psms'

        if not data_name in self.__data_cache.keys():
            ## Get Quan columns
            inspector = inspect(self.__engine)

            quan_value_pattern = re.compile(r'^QuanValue\w+$')
            quan_cols = []
            for column in inspector.get_columns('TargetPsms'):
                if re.match(r'^QuanValue\w+$', column['name']):
                    quan_cols.append(column['name'])

            sqlStr = """
                        SELECT
                        Sequence,
                        ModifiedSequence,
                        Modifications,
                        ParentProteinAccessions,
                        ParentProteinDescriptions,
                        SpectrumFileName,
                        QuanChannel,
                        TargetPeptideGroupsPeptideGroupID,
                        {0}
                        FROM TargetPsms t1, TargetPeptideGroupsTargetPsms t2 
                        WHERE t1.PeptideID = t2.TargetPsmsPeptideID AND QuanChannel IS NOT NULL
                    """.format(','.join(quan_cols))

            self.__data_cache[data_name] = self.__read_data_frame(sqlStr)

        return self.__data_cache[data_name]

    def get_target_peptides(self):
        data_name = 'target_peptides'

        if not data_name in self.__data_cache.keys():
            sqlStr = """
                        SELECT
                        PeptideGroupID, 
                        Checked,
                        Confidence,
                        ExcludedBy,
                        Sequence,
                        Modifications_all_positions,
                        Modifications_best_positions,
                        Contaminant,
                        QvalityPEP,
                        Qvalityqvalue,
                        ParentProteinGroupCount,
                        ParentProteinCount,
                        PsmCount,
                        MasterProteinAccessions,
                        MissedCleavages,
                        TheoreticalMass,
                        QuanInfo,
                        IonsScoreMascot,
                        ConfidenceMascot,
                        PercolatorqValueMascot,
                        PercolatorPEPMascot,
                        AbundanceRatios,
                        Abundances
                        FROM TargetPeptideGroups tpg
                        WHERE AbundanceRatios IS NOT NULL;
                    """

            df = self.__read_data_frame(sqlStr)

            df['local_mod_locs'] = df['Modifications_best_positions'].apply(self.__get_local_mod_locations)
            df['global_mod_locs'] = df.apply(self.__get_global_mod_locations, axis=1)

            df['files'] = df['PeptideGroupID'].apply(self.__get_found_raw_files)

            df['AbundanceRatios'] = df['AbundanceRatios'].apply(self.__extract_values, n=1, dataType=DataType.Float)
            df['Abundances'] = df['Abundances'].apply(self.__extract_values,
                                                      n=self.__num_quant_channels,
                                                      dataType=DataType.Float)

            df['Log2Ratio'] = df['AbundanceRatios'].apply(np.log2)

            columns = ['Sequence',
                       'Modifications_all_positions', 'Modifications_best_positions',
                       'local_mod_locs', 'global_mod_locs',
                       'QvalityPEP', 'Qvalityqvalue', 'ParentProteinGroupCount',
                       'ParentProteinCount', 'PsmCount', 'MasterProteinAccessions',
                       'Checked', 'Confidence', 'ExcludedBy', 'Contaminant',
                       'MissedCleavages', 'TheoreticalMass', 'QuanInfo', 'IonsScoreMascot',
                       'ConfidenceMascot', 'PercolatorqValueMascot', 'PercolatorPEPMascot',
                       'AbundanceRatios', 'Log2Ratio', 'Abundances',
                       'files']

            self.__data_cache[data_name] = df[columns]

        return self.__data_cache[data_name]
