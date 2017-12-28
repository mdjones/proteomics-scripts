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
import xml.etree.ElementTree as ET


class CDataType(Enum):
    #https: // www.tutorialspoint.com / cprogramming / c_data_types.htm
    Double = 1
    Long = 2


class PDReader:
    __peptide_group_mod_pattern = re.compile(r'\[(\w\d+)\]')
    __peptide_group_isotop_mod_pattern = re.compile(r'\d×IsoTOP[\w ]*\[(\w\d+)\]', re.IGNORECASE)
    url = None

    def __init__(self, pd_result_file=None,
                 include_non_quant=False,
                 pd_version='2.1'):
        """

        :param pd_result_file:
        :param include_non_quant:
        :param pd_version: - could get this from  SchemaInfo table in PDResults
        """

        self.url = 'sqlite:///{0}'.format(pd_result_file)
        self.__engine = create_engine(self.url)
        self.__include_non_quant = include_non_quant
        self.__pd_version = pd_version
        self.__data_cache = {}

    def __extract_values(self, binary_data):
        """
        Adapted from Tomas Reijtar Unicorn project.

        Number of values and type of data inferred from binary_data
        See: https://www.tutorialspoint.com/cprogramming/c_data_types.htm
        The binary data length will always be 1+length of data type. The first
        byte represents the Byte Order, Size, and Alignment
        (See: https://docs.python.org/2/library/struct.html)

        So far the types we support are:
        Long - 4 bytes
        Double - 8 bytes

        Parameters
        ----------
        binary_data
                    C Struct represented as a python string (Ex:  b'\xe7rz\xec\xf8\xc2\xe2?\x01')

        Returns
                A single value or an array depending on the number of values in binary_data
        -------

        """

        if binary_data:
            result = []
            if (len(binary_data) % 9) == 0:
                dataType = CDataType.Double
            elif (len(binary_data) % 5) == 0:
                dataType = CDataType.Long
            else:
                raise ValueError('Cannot determine DataType for binary {} of size {}'.format(binary_data, len(binary_data)))

            if dataType == CDataType.Double:
                n = int(len(binary_data)/9)
                for i in range(n):
                    sub = binary_data[9 * i:9 * i + 8]
                    result.append(struct.unpack("d", sub)[0])
            elif dataType == CDataType.Long:
                n = int(len(binary_data)/5)
                for i in range(n):
                    result.append(struct.unpack("i", binary_data[5 * i:5 * i + 4]))
            else:
                raise ValueError('Unsupported DataType {}'.format(dataType))


        # take care of the missing values
        else:
            result = [np.nan]

        if len(result) == 1:
            result = result[0]

        return result

    def __get_found_raw_files(self, targetPeptideGroupsPeptideGroupID):
        df = self.__get_peptidegroupid_to_spectrumfile()
        df = df[df['TargetPeptideGroupsPeptideGroupID'] == targetPeptideGroupsPeptideGroupID]

        values = df['SpectrumFileName'].values

        return values

    def __quantile_normalized_ratio(self, abundances, max_ratio=100):
        # Taken from
        # https://stackoverflow.com/questions/37935920/quantile-normalization-on-pandas-dataframe

        df = abundances.apply(pd.Series)
        assert df.columns.size == 2, 'Can only take ratios of dataframes with 2 columns'

        rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
        df = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()

        min_ratio = 1/max_ratio
        ratio = (df[0]/df[1]).replace(np.Infinity, max_ratio).replace(0, min_ratio)

        return ratio

    def __get_local_mod_locations(self, modifications, patern=None):
        mod_locs = re.findall(patern, modifications)
        return mod_locs

    def __get_global_mod_locations(self, row, mod_col=None):
        local_mods = row[mod_col]
        local_mods = [re.match(r'(\w)(\d+)', mod).groups() for mod in local_mods]
        peptideGroupID = row['PeptideGroupID']
        peptide_sequence = row['Sequence']
        accessions = row['MasterProteinAccessions']
        if accessions:
            accessions = accessions.split('; ')
        else:
            accessions = []

        df = self.__get_target_proteins()
        df = df[df['TargetPeptideGroupsPeptideGroupID'] == peptideGroupID]

        peptide_starts = []
        for accession in accessions:
            proteinSequence = df[df['Accession'] == accession]['Sequence'].values
            assert len(proteinSequence) <= 1, 'More then one sequence returned for {} -- {}'.format(accession, proteinSequence)
            if len(proteinSequence) == 1:
                peptide_starts.append((proteinSequence[0].find(peptide_sequence)))

        global_positions = []
        for peptide_start in peptide_starts:
            pos_strings = []

            for local_mod in local_mods:
                local_position = int(local_mod[1])
                global_pos = np.nan if peptide_start == -1 else (local_position + peptide_start)
                pos_strings.append('{0}{1}'.format(local_mod[0],global_pos))

            global_positions.append(','.join(pos_strings))

        return global_positions

    def __create_protein_isotop_locs(self, row):
        locs = row['GLOBAL_ISOTOP_LOCS']
        accessions = row['MasterProteinAccessions']
        if accessions:
            accessions = accessions.split('; ')
        else:
            accessions = []

        return ','.join([a + '_' + b for a, b in zip(accessions, locs)])

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

    def __get_peptidegroupid_to_spectrumfile(self):
        data_name = 'peptide_to_file'
        if not data_name in self.__data_cache.keys():
            quan_channel_filter = '' if self.__include_non_quant else ' AND QuanChannel IS NOT NULL'

            sqlStr = """
                        SELECT
                        SpectrumFileName,TargetPeptideGroupsPeptideGroupID
                        FROM TargetPsms t1,TargetPeptideGroupsTargetPsms t2
                        WHERE t1.PeptideID = t2.TargetPsmsPeptideID {0}
                        GROUP BY SpectrumFileName,TargetPeptideGroupsPeptideGroupID
                    """.format(quan_channel_filter)

            self.__data_cache[data_name] = self.__read_data_frame(sqlStr)

        return self.__data_cache[data_name]

    def get_target_peptides(self, include_additional_data=True):
        targetPeptideGroupsTableDF = self.get_target_peptide_groups_table()
        if include_additional_data:
            targetPeptideGroupsTableDF = targetPeptideGroupsTableDF.drop(labels=['AbundanceRatios', 'Abundances'], axis=1)

            additionalTargetPepDF = self.get_additional_target_peptide_data()
            df = pd.merge(targetPeptideGroupsTableDF,
                          additionalTargetPepDF,
                          left_index=True,
                          right_index=True,
                          indicator=True)

            assert df[df._merge != 'both'].empty, 'Invalid merge {}'.format(df[df['_merge'] != 'both'])

            df = df.drop(labels=['_merge'], axis=1)
            return df
        else:
            return targetPeptideGroupsTableDF

    def get_analysis_definition(self):
        data_set_name = 'analysis_definition'
        if not data_set_name in self.__data_cache.keys():
            sqlStr = """
                                Select AnalysisDefinitionXML FROM AnalysisDefinition"""

            self.__data_cache[data_set_name] = self.__read_data_frame(sqlStr).iloc[0,0]

        return self.__data_cache[data_set_name]

    def get_target_peptide_groups_table(self):
        quan_channel_filter = '' if self.__include_non_quant else ' WHERE AbundanceRatios IS NOT NULL'

        data_set_name = 'target_peptides'
        if not data_set_name in self.__data_cache.keys():
            sqlStr = """
                        SELECT
                        PeptideGroupID,
                        Checked,
                        Confidence,
                        ExcludedBy,
                        Sequence,
                        Modifications_all_positions,
                        Modifications_best_positions,
                        QvalityPEP,
                        Qvalityqvalue,
                        ParentProteinGroupCount,
                        ParentProteinCount,
                        PsmCount,
                        MasterProteinAccessions,
                        MissedCleavages,
                        TheoreticalMass,
                        QuanInfo,
                        AbundanceRatios,
                        Abundances
                        FROM TargetPeptideGroups tpg
                        {0}
                """.format(quan_channel_filter)

            self.__data_cache[data_set_name] = self.__read_data_frame(sqlStr)

        return self.__data_cache[data_set_name]

    def get_additional_target_peptide_data(self):
        data_set_name = 'target_peptides_plus'
        if not data_set_name in self.__data_cache.keys():
            targetPeptideDF = self.get_target_peptide_groups_table()
            df = pd.DataFrame(index=targetPeptideDF.index)

            df['FILES'] = targetPeptideDF['PeptideGroupID'].apply(self.__get_found_raw_files)
            df['ABUNDANCE_RATIOS'] = targetPeptideDF['AbundanceRatios'].apply(self.__extract_values)
            df['ABUNDANCE_LOG2_RATIO'] = df['ABUNDANCE_RATIOS'].apply(np.log2)

            df['ABUNDANCES'] = targetPeptideDF['Abundances'].apply(self.__extract_values)
            named_abundances = self.__get_named_abundances(df['ABUNDANCES'])
            for col in named_abundances.columns:
                df[col] = named_abundances[col]

            df['QUANTILE_NORM_ABUNDANCE_RATIO']  = self.__quantile_normalized_ratio(df['ABUNDANCES'])
            df['QUANTILE_NORM_ABUNDANCE_LOG2_RATIO'] = df['QUANTILE_NORM_ABUNDANCE_RATIO'].apply(np.log2)

            df['Sequence'] = targetPeptideDF['Sequence']
            df['PeptideGroupID'] = targetPeptideDF['PeptideGroupID']
            df['MasterProteinAccessions'] = targetPeptideDF['MasterProteinAccessions']

            df['MOD_LOCS'] = targetPeptideDF['Modifications_best_positions'].apply(
                lambda x : self.__get_local_mod_locations(x, patern=self.__peptide_group_mod_pattern))
            df['GLOBAL_MOD_LOCS'] = df.apply(
                lambda x : self.__get_global_mod_locations(x, mod_col='MOD_LOCS'), axis=1)

            df['ISOTOP_MOD_LOCS'] = targetPeptideDF['Modifications_best_positions'].apply(
                lambda x: self.__get_local_mod_locations(x, patern=self.__peptide_group_isotop_mod_pattern))
            df['GLOBAL_ISOTOP_LOCS'] = df.apply(
                lambda x : self.__get_global_mod_locations(x, mod_col='ISOTOP_MOD_LOCS'), axis=1)

            df['protein_isotop_locs'] = df.apply(self.__create_protein_isotop_locs, axis=1)
            df = df.drop(labels=['Sequence', 'PeptideGroupID', 'MasterProteinAccessions'], axis=1)


            self.__data_cache[data_set_name] = df

        return self.__data_cache[data_set_name]

    def __get_named_abundances(self, abundances_df):
        xml_def = self.get_analysis_definition()
        root = ET.fromstring(xml_def)

        channels = []
        channel_num = 1
        for child in root.findall('./StudyDefinition/QuanMethods/QuanMethod/QuanChannels/QuanChannel'):
            child_name = '{}_{}'.format(child.attrib['Name'], child.attrib['Position']).replace(' ', '_').upper()
            channels.append(child_name)

        df = abundances_df.apply(pd.Series)
        df.columns = channels

        return df
