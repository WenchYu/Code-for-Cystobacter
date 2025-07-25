
'''
MSanalyst annotating non-targeted metabolomics data includes the following steps:
1. MS1 match conducted by 'ms1_match()'
2. 'ISDB_MS2_match()' and 'EDB_MS2_match()' will perform MS2 comparison based on the MS1 match results
3. Query features are first clustered and curated by 'self_clustering()'
4. The network will be annotated by filtered MS2 comparison results in 'molecular_generation()'
'''

import os
import ast,time,json,heapq,spectral_entropy
import pandas as pd
import numpy as np
import networkx as nx
import spectrum_utils.spectrum as sus
from tqdm import tqdm, trange
from joblib import Parallel, delayed
from spectral_entropy import similarity
from my_packages import functions
from my_packages.peaktools import neutral_loss,modified_cosine,find_match_peaks_efficient,convert_to_peaks


class FlashPrecursorSearch:
    def __init__(
            self,
            max_ms2_tolerance_in_da=0.02,
            mz_index_step=0.0001,
            adduct='m+h') -> None:
        """
        Initialize the EntropySearch class.

        :param max_ms2_tolerance_in_da: The maximum MS2 tolerance used when searching the MS/MS spectra, in Dalton. Default is 0.024.
        :param mz_index_step:   The step size of the m/z index, in Dalton. Default is 0.0001.
                                The smaller the step size, the faster the search, but the larger the index size and longer the index building time.
        """
        self.mz_index_step = mz_index_step
        self.max_ms2_tolerance_in_da = max_ms2_tolerance_in_da
        self.adduct = adduct

        self.index = []  # record ions_mz

        self.index_names = [
            "all_ions_mz_idx_start",
            "all_ions_mz",
        ]
        self.index_dtypes = {
            "all_ions_mz_idx_start": np.int64,
            "all_ions_mz": np.float32,
        }

    def _generate_index_from_peak_data(self, dataframe, max_indexed_mz):
        # Sort with precursor m/z and pre-sort

        # Record the precursor m/z
        all_ions_mz = dataframe[self.adduct]

        # Build index for fast access to the ion's m/z.
        max_mz = max_indexed_mz
        search_array = np.arange(0.0, max_mz, self.mz_index_step)
        all_ions_mz_idx_start = np.searchsorted(all_ions_mz, search_array, side="left").astype(np.int64)

        ############## Step 4: Save the index. ##############
        index = [
            all_ions_mz_idx_start,
            all_ions_mz,
        ]
        return index

    def build_index(self, dataframe, max_indexed_mz: float = 1500.00005):
        """
        Build the index for the MS/MS spectra library.

        The spectra provided to this function should be a dictionary in the format of {"precursor_mz": precursor_mz, "peaks": peaks}.
        The precursor_mz is the precursor m/z value of the MS/MS spectrum;
        The peaks is a numpy array which has been processed by the function "clean_spectrum".

        :param all_spectra_list:    A list of dictionaries in the format of {"precursor_mz": precursor_mz, "peaks": peaks},
                                    the spectra in the list need to be sorted by the precursor m/z.
        :param max_indexed_mz: The maximum m/z value that will be indexed. Default is 1500.00005.
        """

        ############## Step 2: Build the index by sort with product ions. ##############
        self.index = self._generate_index_from_peak_data(dataframe, max_indexed_mz)
        return self.index

    def _find_location_from_array_with_index(self, wanted_mz, mz_array, mz_idx_start_array, side,
                                             index_number_in_one_da):
        mz_min_int = (np.floor(wanted_mz * index_number_in_one_da)).astype(int)
        mz_max_int = mz_min_int + 1

        if mz_min_int >= len(mz_idx_start_array):
            mz_idx_search_start = mz_idx_start_array[-1]
        else:
            mz_idx_search_start = mz_idx_start_array[mz_min_int].astype(int)

        if mz_max_int >= len(mz_idx_start_array):
            mz_idx_search_end = len(mz_array)
        else:
            mz_idx_search_end = mz_idx_start_array[mz_max_int].astype(int) + 1

        return mz_idx_search_start + np.searchsorted(mz_array[mz_idx_search_start:mz_idx_search_end], wanted_mz,
                                                     side=side)

    def _get_indeces(self, mz_query):
        index_number_in_one_da = int(1/self.mz_index_step)
        all_ions_mz = self.index[1]
        all_ions_mz_idx_start = self.index[0]
        product_mz_idx_min = self._find_location_from_array_with_index(
            mz_query - self.max_ms2_tolerance_in_da, all_ions_mz, all_ions_mz_idx_start, "left",
            index_number_in_one_da)

        product_mz_idx_max = self._find_location_from_array_with_index(
            mz_query + self.max_ms2_tolerance_in_da, all_ions_mz, all_ions_mz_idx_start, "right",
            index_number_in_one_da)

        return product_mz_idx_min, product_mz_idx_max

def match_mz(quant_df_row, msdb_df, mz_column='row m/z',isdb_ms1_match_threshld = 5):
    '''
    Single MS1 match against the in-silico MS1 library

    :param quant_df_row: Rows in quantification table(xxx_quant.csv) generated by MZmine
    :param msdb_df: The default database is (default: '../msdb/isdbMS1.csv')
    :param mz_column: Column name of queried MS1 (default: 'row m/z')
    :param np_ms1_match_threshold: Relative error (ppm)
    :return:tuple([...],[...],[...])
    '''
    hits_id = []
    hits_smiles = []
    for j in range(len(msdb_df.index)):
        ppm_H = functions.calculate_ppm(quant_df_row[mz_column], msdb_df.loc[j,'m+h'])
        ppm_Na = functions.calculate_ppm(
                quant_df_row[mz_column], msdb_df.loc[j,'m+na'])

        if ppm_H < isdb_ms1_match_threshld or ppm_Na < isdb_ms1_match_threshld:
            hits_id.append(msdb_df.id[j])
            hits_smiles.append(msdb_df.smiles[j])
    if not hits_id or not hits_smiles:
        hits_id.append(None)
        hits_smiles.append(None)
    return hits_id, hits_smiles

def match_edb_mz(quant_df_row,edb_df,mz_column='row m/z',edb_ms1_match_threshold = 5):
    '''
    Single MS1 match against experimental MS1 library
    '''
    hits_id = []
    hits_smiles = []
    for j in range(len(edb_df.index)):
        ppm = functions.calculate_ppm(quant_df_row[mz_column], edb_df['pepmass'][j])
        if ppm < edb_ms1_match_threshold:
            hits_id.append(edb_df.id[j])
            hits_smiles.append(str(edb_df.smiles[j]))
    if not hits_id or not hits_smiles:
        hits_id.append(None)
        hits_smiles.append(None)
    return hits_id, hits_smiles

def ms1_match(args,queryDF = None):
    '''
    MS1 match against the entire MS1 library (experimental and in-silico)
    :param args: args.pepmass_match_tolerance,
    :return: Two MS1 match result files, prefixed with 'IS_MS1match_' and 'E_MS1match_'
    '''

    isdb_ppm = args.pepmass_match_tolerance # Setting allowed mass error
    edb_ppm = args.pepmass_match_tolerance

    isdb_df = pd.read_csv(args.isms1_file,low_memory=False) # Loading MS1 library
    edb_df = functions.df_preprocess(args.edbms1_file)

    if queryDF is None:
        query = args.quant_file # Query dataframe from quantification table
        quant_df = functions.df_preprocess(query) # loading quantification table
        quant_name = os.path.splitext(os.path.basename(args.quant_file))[0]
        result_dir = os.path.join(args.output, f'{quant_name}_result')
        os.makedirs(result_dir, exist_ok=True)
        basename = os.path.basename(args.quant_file)
    else:
        quant_df= queryDF # Query dataframe from single pepmass (used in '../ms1search.py')
        query_mz = str(queryDF.iloc[0,1])
        result_dir = os.path.join(args.output,f'{query_mz}')
        os.makedirs(result_dir, exist_ok=True)
        basename = f'{query_mz}.csv'

    # MS1 match against the entire library (experimental and in-silico)
    n_jobs = args.cpus
    isdb_results = Parallel(n_jobs=n_jobs)(
        delayed(match_mz)(quant_df_row, isdb_df, isdb_ms1_match_threshld=isdb_ppm) for quant_df_row in
        tqdm(quant_df.to_dict('records')))
    edb_results = Parallel(n_jobs=n_jobs)(
        delayed(match_edb_mz)(quant_df_row, edb_df, edb_ms1_match_threshold=edb_ppm) for quant_df_row in
        tqdm(quant_df.to_dict('records')))

    # Output
    isdb_match_rows = []
    edb_match_rows = []
    for i, (hits_id, hits_smiles) in enumerate(isdb_results):
        for j in range(len(hits_id)):
            isdb_match_row = {'row ID': quant_df.at[i, 'row ID'], 'row m/z': quant_df.at[i, 'row m/z'],
                            'match_id': hits_id[j], 'match_smiles': hits_smiles[j]}
            isdb_match_rows.append(isdb_match_row)
    isdb_match_df = pd.DataFrame(isdb_match_rows)
    for i, (hits_id, hits_smiles) in enumerate(edb_results):
        for j in range(len(hits_id)):
            edb_match_row = {'row ID': quant_df.at[i, 'row ID'], 'row m/z': quant_df.at[i, 'row m/z'],
                              'match_id': hits_id[j], 'match_smiles': hits_smiles[j]}
            edb_match_rows.append(edb_match_row)
    edb_match_df = pd.DataFrame(edb_match_rows)

    isdb_result_path = os.path.join(result_dir, f'IS_MS1match_{os.path.basename(basename)}')
    edb_result_path = os.path.join(result_dir, f'E_MS1match_{os.path.basename(basename)}')

    isdb_match_df.to_csv(isdb_result_path, index=False)
    edb_match_df.to_csv(edb_result_path, index=False)

    quant_df['isdbms1_id'] = np.nan
    quant_df['isdbms1_smiles'] = np.nan
    quant_df['edbms1_id'] = np.nan
    quant_df['edbms1_smiles'] = np.nan
    for i, (hits_id, hits_smiles) in enumerate(isdb_results):
        if hits_id is not None and all(isinstance(x, str) for x in hits_id):
            quant_df.at[i, 'isdbms1_id'] = ';'.join([x or '' for x in hits_id])
            quant_df.at[i, 'isdbms1_smiles'] = ';'.join(hits_smiles)
    for i, (hits_id, hits_smiles) in enumerate(edb_results):
        if hits_id is not None and all(isinstance(x, str) for x in hits_id):
            quant_df.at[i, 'edbms1_id'] = ';'.join([x or '' for x in hits_id])
            quant_df.at[i, 'edbms1_smiles'] = ';'.join(hits_smiles)

    ms1_result_path = os.path.join(result_dir, f'MS1match_{basename}')
    quant_df.to_csv(ms1_result_path, index=False)
    print('MS1 matching finished!')

def ISDB_MS2_match(args,queryMGF=None):
    ''' MS2 match against in-silico MS2 library '''
    if queryMGF is None:
        mgf_file = args.mgf_file # Query MS2 from '.mgf' file
        quant_name = os.path.splitext(os.path.basename(args.quant_file))[0]
        result_dir = os.path.join(args.output, f'{quant_name}_result')
        os.makedirs(result_dir, exist_ok=True)
        basename = os.path.basename(args.quant_file)
        exp_info = functions.mgf_process(mgf_file)
    else:
        spectra_info = queryMGF # Query MS2 from direct input (used in '../ms2search.py')
        query_mz = spectra_info.loc[0,'pepmass']
        result_dir = os.path.join(args.output, f'{query_mz}')
        os.makedirs(result_dir, exist_ok=True)
        basename = f'{query_mz}.csv'
        exp_info = queryMGF


    with open(args.isms2_file) as f:
        isdb_info = json.load(f) # Loading in-silico MS2 library

    # Loading in-silico MS1 match result (csv file generated by ms1_match())
    is_result_path = os.path.join(result_dir, f'IS_MS1match_{basename}')
    is_ms1_match_df = functions.df_preprocess(is_result_path)

    # 0, 1, 2 represent the in-silico MS2 predicted by CFM-ID at three collision energy (10, 20 and 40 eV)
    is_ms1_match_df['mps0'] = np.nan  # mps: matched peaks
    is_ms1_match_df['pp0'] = np.nan # pp: peak percentage
    is_ms1_match_df['pair_similarity0'] = np.nan # MS2 spectral similarity
    is_ms1_match_df['mps1'] = np.nan
    is_ms1_match_df['pp1'] = np.nan
    is_ms1_match_df['pair_similarity1'] = np.nan
    is_ms1_match_df['mps2'] = np.nan
    is_ms1_match_df['pp2'] = np.nan
    is_ms1_match_df['pair_similarity2'] = np.nan

    for i in trange(len(is_ms1_match_df)):
        row_id = str(is_ms1_match_df.loc[i,'row ID']) # convert int to str
        match_id = str(is_ms1_match_df.loc[i,'match_id'])
        if match_id != 'nan':
            try:  # Some features in xxx_quant.csv have no MS2
                exp_pm = float(exp_info[exp_info['id'] == row_id].pepmass.iloc[0]) # pepmass of query feature
                exp_charge = int(exp_info[exp_info['id'] == row_id].charge.iloc[0]) # charge of query feature
                exp_ms2 = np.asarray(exp_info[exp_info['id'] == row_id].ms2.iloc[0]) # MS2 of query feature
                exp_ms2 = spectral_entropy.clean_spectrum(exp_ms2, max_mz=exp_pm + 0.01) # MS2 spectrum clean by normalizing and removing signals with intensity less than 1% of the base peak
                exp_mz = np.array(exp_ms2[:, 0], dtype=np.float64)
                exp_intensity = np.array(exp_ms2[:, 1], dtype=np.float64)
                exp_spectrum = sus.MsmsSpectrum(identifier=row_id, precursor_mz=exp_pm + 0.01, precursor_charge=exp_charge,
                                                mz=exp_mz,intensity=exp_intensity) # Spectrum for MS2 comparison using modified_cosine and neutral_loss

                is_pm = exp_pm
                is_smile = isdb_info[match_id]['smiles']
                # 0, 1, 2 represent the in-silico MS2 predicted by CFM-ID at three collision energy (10, 20 and 40 eV)
                e0_ms2 = np.asarray(ast.literal_eval(isdb_info[match_id]['energy0_ms2']))
                e0_mz = np.array(e0_ms2[:,0],dtype=np.float64)
                e0_intensity = np.array(e0_ms2[:, 0], dtype=np.float64)
                e0_spectrum = sus.MsmsSpectrum(identifier=f'e0_{match_id}', precursor_mz=is_pm + 0.01, precursor_charge=1,
                                                mz=e0_mz,intensity=e0_intensity)
                e1_ms2 = np.asarray(ast.literal_eval(isdb_info[match_id]['energy1_ms2']))
                e1_mz = np.array(e1_ms2[:, 0], dtype=np.float64)
                e1_intensity = np.array(e1_ms2[:, 0], dtype=np.float64)
                e1_spectrum = sus.MsmsSpectrum(identifier=f'e1_{match_id}', precursor_mz=is_pm + 0.01, precursor_charge=1,
                                               mz=e1_mz,intensity=e1_intensity)

                e2_ms2 = np.asarray(ast.literal_eval(isdb_info[match_id]['energy2_ms2']))
                e2_mz = np.array(e2_ms2[:, 0], dtype=np.float64)
                e2_intensity = np.array(e2_ms2[:, 0], dtype=np.float64)
                e2_spectrum = sus.MsmsSpectrum(identifier=f'e2_{match_id}', precursor_mz=is_pm + 0.01, precursor_charge=1,
                                               mz=e2_mz,
                                               intensity=e2_intensity)

                shift = abs(is_pm-exp_pm) # Allowed mass shift in MS2 comparison using modified_cosine
                exp_peaks = len(exp_ms2) # Number of fragments of query feature
                sim0, sim1, sim2 = 0.0, 0.0, 0.0
                mps0, mps1, mps2 = 0, 0, 0
                pp0, pp1, pp2 = 0.0, 0.0, 0.0
                if args.library_matching_method == 'modified_cosine':
                    try:
                        result0 = modified_cosine(exp_spectrum, e0_spectrum, fragment_mz_tolerance=0.05)
                        sim0 = result0.score
                        mps0 = result0.matches
                        pp0 = mps0/exp_peaks
                        result1 = modified_cosine(exp_spectrum, e1_spectrum, fragment_mz_tolerance=0.05)
                        sim1 = result1.score
                        mps1 = result1.matches
                        pp1 = mps1 / exp_peaks
                        result2 = modified_cosine(exp_spectrum, e2_spectrum, fragment_mz_tolerance=0.05)
                        sim2 = result2.score
                        mps2 = result2.matches
                        pp2 = mps2 / exp_peaks
                    except:
                        pass

                elif args.library_matching_method == 'neutral_loss':
                    try:
                        result0 = neutral_loss(exp_spectrum, e0_spectrum, fragment_mz_tolerance=0.05)
                        sim0 = result0.score
                        mps0 = result0.matches
                        pp0 = mps0 / exp_peaks
                        result1 = neutral_loss(exp_spectrum, e1_spectrum, fragment_mz_tolerance=0.05)
                        sim1 = result1.score
                        mps1 = result1.matches
                        pp1 = mps1 / exp_peaks
                        result2 = neutral_loss(exp_spectrum, e2_spectrum, fragment_mz_tolerance=0.05)
                        sim2 = result2.score
                        mps2 = result2.matches
                        pp2 = mps2 / exp_peaks
                    except:
                        pass
                else:
                    mps0 = len(find_match_peaks_efficient(convert_to_peaks(exp_ms2)
                                                          , convert_to_peaks(e0_ms2), shift, 0.05))
                    sim0 = similarity(exp_ms2, e0_ms2, method=args.library_matching_method, ms2_da=0.05)
                    pp0 = mps0 / exp_peaks
                    mps1 = len(find_match_peaks_efficient(convert_to_peaks(exp_ms2)
                                                          , convert_to_peaks(e1_ms2), shift, 0.05))
                    sim1 = similarity(exp_ms2, e1_ms2, method=args.library_matching_method, ms2_da=0.05)
                    pp1 = mps1 / exp_peaks
                    mps2 = len(find_match_peaks_efficient(convert_to_peaks(exp_ms2)
                                                          , convert_to_peaks(e2_ms2), shift, 0.05))
                    sim2 = similarity(exp_ms2, e2_ms2, method=args.library_matching_method, ms2_da=0.05)
                    pp2 = mps2 / exp_peaks

                # Output
                is_ms1_match_df.loc[i, 'pair_similarity0'] = sim0
                is_ms1_match_df.loc[i, 'mps0'] = mps0
                is_ms1_match_df.loc[i, 'pp0'] = pp0
                is_ms1_match_df.loc[i, 'pair_similarity1'] = sim1
                is_ms1_match_df.loc[i, 'mps1'] = mps1
                is_ms1_match_df.loc[i, 'pp1'] = pp1
                is_ms1_match_df.loc[i, 'pair_similarity2'] = sim2
                is_ms1_match_df.loc[i, 'mps2'] = mps2
                is_ms1_match_df.loc[i, 'pp2'] = pp2
                is_ms2_path = os.path.join(result_dir, row_id, f'{match_id}.mgf')
                with open(is_ms2_path, 'w') as f:
                    f.write('BEGIN IONS\n')
                    f.write(f'ID={match_id}\n')
                    f.write(f'PEPMASS={is_pm}\n')
                    f.write(f'SMILES={is_smile}\n')
                    f.write('ENERGY\n')
                    for item in e0_ms2:
                        f.write("%s %s\n" % (item[0], item[1]))
                    f.write('ENERGY1\n')
                    for item in e1_ms2:
                        f.write("%s %s\n" % (item[0], item[1]))
                    f.write('ENERGY2\n')
                    for item in e2_ms2:
                        f.write("%s %s\n" % (item[0], item[1]))
                    f.write('END IONS\n')
            except:
                pass
    is_ms1_match_df.to_csv(is_result_path)

def EDB_MS2_match(args,queryMGF=None):
    ''' MS2 match against experimental MS2 library '''
    if queryMGF is None:
        exp_info = functions.mgf_process(args.mgf_file)  # Query MS2 from '.mgf' file
        quant_name = os.path.splitext(os.path.basename(args.quant_file))[0]
        result_dir = os.path.join(args.output, f'{quant_name}_result')
        os.makedirs(result_dir, exist_ok=True)
        basename = os.path.basename(args.quant_file)
    else:
        spectra_info = queryMGF # Query MS2 from direct input (used in '../ms2search.py')
        query_mz = spectra_info.loc[0,'pepmass']
        result_dir = os.path.join(args.output, f'{query_mz}')
        os.makedirs(result_dir, exist_ok=True)
        basename = f'{query_mz}.csv'
        exp_info = queryMGF

    with open (args.edbms2_file,'r') as f:
        edbms2_info=json.load(f) # Loading experimental MS2 library
    # Loading experimental MS1 match result (csv file generated by ms1_match())
    edb_result_path = os.path.join(result_dir, f'E_MS1match_{basename}')
    edb_ms1_df = functions.df_preprocess(edb_result_path)

    edb_ms1_df['mps'] = np.nan
    edb_ms1_df['pair_similarity'] = np.nan
    edb_ms1_df['pp'] = np.nan
    for i in trange(len(edb_ms1_df)):
        row_id = str(edb_ms1_df.loc[i,'row ID'])
        match_id = str(edb_ms1_df.loc[i,'match_id'])

        if match_id != 'nan':
            try:  # some features have no MS2
                exp_pm = float(exp_info[exp_info['id'] == row_id].pepmass.iloc[0]) # pepmass of query feature
                exp_ms2 = exp_info[exp_info['id'] == row_id].ms2.iloc[0] # ms2 of query feature
                exp_ms2 = spectral_entropy.clean_spectrum(exp_ms2, max_mz=exp_pm+0.01) # MS2 spectrum clean by normalizing and removing signals with intensity less than 1% of the base peak
                exp_charge =  int(exp_info[exp_info['id'] == row_id].charge.iloc[0]) # charge of query feature
                exp_mz = np.array(exp_ms2[:, 0], dtype=np.float64)
                exp_intensty = np.array(exp_ms2[:, 1], dtype=np.float64)
                exp_spectrum = sus.MsmsSpectrum(identifier=row_id, precursor_mz=exp_pm+0.01, precursor_charge=exp_charge, mz=exp_mz,
                                             intensity=exp_intensty) # Spectrum for MS2 comparison using modified_cosine and neutral_loss

                edb_pm = float(edbms2_info[match_id]['pepmass'])
                edb_smiles = edbms2_info[match_id]['smiles']
                edb_ms2 = np.asarray(ast.literal_eval(edbms2_info[match_id]['ms2']))
                edb_ms2 = spectral_entropy.clean_spectrum(edb_ms2, max_mz=edb_pm+0.01)
                edb_charge = int(edbms2_info[match_id]['charge'])
                edb_mz = np.array(edb_ms2[:, 0], dtype=np.float64)
                edb_intensty = np.array(edb_ms2[:, 1], dtype=np.float64)
                edb_spectrum = sus.MsmsSpectrum(identifier=match_id, precursor_mz=edb_pm+0.01, precursor_charge=edb_charge, mz=edb_mz,
                                             intensity=edb_intensty)


                shift = abs(exp_pm-edb_pm) # Allowed mass shift in MS2 comparison using modified_cosine
                exp_peaks = len(exp_ms2) # Number of fragments of query feature
                sim, mps, pp = 0.0, 0, 0.0
                if args.library_matching_method == 'modified_cosine':
                    try:
                        result = modified_cosine(exp_spectrum, edb_spectrum, fragment_mz_tolerance=0.05)
                        sim = result.score
                        mps = result.matches
                        pp = mps/exp_peaks
                    except:
                        pass
                elif args.library_matching_method == 'neutral_loss':
                    try:
                        result = neutral_loss(exp_spectrum, edb_spectrum, fragment_mz_tolerance=0.05)
                        sim = result.score
                        mps = result.matches
                        pp = mps / exp_peaks
                    except:
                        pass
                else:
                    mps = len(find_match_peaks_efficient(convert_to_peaks(exp_ms2)
                                                          , convert_to_peaks(edb_ms2), shift, 0.05))
                    sim = similarity(exp_ms2, edb_ms2, method=args.library_matching_method, ms2_da=0.05)
                    pp = mps/exp_peaks

                # Output
                edb_ms1_df.loc[i, 'pair_similarity'] = sim
                edb_ms1_df.loc[i, 'mps'] = mps
                edb_ms1_df.loc[i,'pp'] = pp
                edb_ms2_path = os.path.join(result_dir, row_id, f'{match_id}.mgf')
                with open(edb_ms2_path, 'w') as f:
                    f.write('BEGIN IONS\n')
                    f.write(f'ID={match_id}\n')
                    f.write(f'PEPMASS={edb_pm}\n')
                    f.write(f'SMILES={edb_smiles}\n')
                    f.write('ENERGY\n')
                    for item in edb_ms2:
                        f.write("%s %s\n" % (item[0], item[1]))
                    f.write('END IONS\n')
            except:
                pass
    edb_ms1_df.to_csv(edb_result_path)
    print('MS2 matching finished!')

def mn_curating(G, topk):
    '''
    Limit the size of the cluster by setting the number of neighbors (topk) allowed for a node

    :param G: Network graph
    :param topk: Maximum neighbors allowed for a node
    :return: An curated G
    '''
    node_ids = list(G.nodes())
    for node_id in node_ids:
        if len(G[node_id]) > topk:
            edges = list(G.edges(node_id)) # List all the edges
            # Keep the topK most similar neighbors of a node
            result = \
                heapq.nlargest(topk,
                               [(data.get('pair_similarity', 0), neighbor) for neighbor, data in G[node_id].items()])
            topk_edges = [t[1] for t in result]
            for edge in edges:
                if edge[1] not in topk_edges:
                    G.remove_edge(edge[0], edge[1])
    return G

def self_clustering(args):
    '''

    :param args: args.output, args.quant_file, args.mgf_file, args.self_clustering_similarity
    :return:
    '''
    parent_folder = f'{args.output}/{os.path.splitext(os.path.basename(args.quant_file))[0]}_result/merge/'
    if not os.path.exists(parent_folder): # Make sure the parent folder exist
        os.makedirs(parent_folder)
    exp_info = functions.mgf_process(args.mgf_file) # Loading query '.mgf' file

    G = nx.MultiGraph()  # Creating undirected graph
    for i, (id1, pm1, charge1, spec1) in exp_info.iterrows():
        pm1 = float(pm1)
        node_attr = {'pepmass': pm1}
        G.add_node(id1, **node_attr)  # add nodes and attributes

    # Self clustering
    for i, (id1, pm1, charge1, spec1) in tqdm(exp_info.iterrows(), total=len(exp_info)):
        pm1 = float(pm1)
        charge1 = int(charge1)
        spec1 = spectral_entropy.clean_spectrum(spec1, max_mz=pm1 - 0.01, noise_removal=0.01)
        mz1 = np.array(spec1[:, 0], dtype=np.float64)
        intensity1 = np.array(spec1[:, 1], dtype=np.float64)
        spectrum1 = sus.MsmsSpectrum(identifier=id1, precursor_mz=pm1, precursor_charge=charge1
                                     , mz=mz1, intensity=intensity1).remove_precursor_peak(0.01, "Da")
        peaks1 = len(spec1)
        G.add_node(id1, **{'num_fragments': peaks1})
        for j, (id2, pm2, charge2, spec2) in exp_info.iloc[:i, ].iterrows():
            pm2 = float(pm2)
            charge2 = int(charge2)
            spec2 = spectral_entropy.clean_spectrum(spec2, max_mz=pm2 - 0.01, noise_removal=0.01)
            mz2 = np.array(spec2[:, 0], dtype=np.float64)
            intensity2 = np.array(spec2[:, 1], dtype=np.float64)
            spectrum2 = sus.MsmsSpectrum(identifier=id2, precursor_mz=pm2, precursor_charge=charge2, mz=mz2,
                                         intensity=intensity2).remove_precursor_peak(0.01, "Da")

            sim, mps, pp = 0.0, 0, 0.0
            if args.self_clustering_method == 'modified_cosine':
                try:
                    result = modified_cosine(spectrum1, spectrum2, fragment_mz_tolerance=0.02)
                    sim = result.score
                    mps = result.matches
                    pp = mps / peaks1
                except:
                    pass
            elif args.self_clustering_method == 'neutral_loss':
                try:
                    result = neutral_loss(spectrum1, spectrum2, fragment_mz_tolerance=0.02)
                    sim = result.score
                    mps = result.matches
                    pp = mps / peaks1
                except:
                    pass
            else:
                try:
                    sim = similarity(spec1, spec2, method=args.self_clustering_method, ms2_da=0.02)
                    result = modified_cosine(spectrum1, spectrum2, fragment_mz_tolerance=0.02)
                    mps = result.matches
                    pp = mps / peaks1
                except:
                    pass
            if sim >= args.self_clustering_similarity \
                    and mps >= args.self_clustering_peaks:
                edge_attr = {'pair_similarity': sim, 'matched_peaks': mps, 'peak_percentage': pp,
                             'edge_type': args.self_clustering_method}
                G.add_edge(id1, id2, **edge_attr)
    # Output
    G = mn_curating(G,args.top_k)
    print('Self clustering finished!')
    MN_file = os.path.join(parent_folder,
                           f'{os.path.splitext(os.path.basename(args.mgf_file))[0]}_{args.self_clustering_method}_{args.self_clustering_similarity}_{args.self_clustering_peaks}.graphml')
    nx.write_graphml(G, MN_file)

def molecular_generation(args):
    '''
    Generating the final molecular network
    A : features without any MS1 match
    B1 : features have MS1 matches in edb or isdb but MS2 match unwell
    B2 : features have MS1 matches in edb or isdb and MS2 match well
    C1 : features have MS1 matches in isdb but MS2 match unwell
    C2 : features have MS1 matches in isdb and MS2 match well
    '''

    parent_folder = f'{args.output}/{os.path.splitext(os.path.basename(args.quant_file))[0]}_result'
    quant_df = functions.df_preprocess(args.quant_file)
    exp_info = functions.mgf_process(args.mgf_file)
    row_ids = [int(x) for x in exp_info['id'].values.tolist()]

    G = nx.MultiGraph()  # Creating undirected graph
    for i, (id1, pm1, charge1, spec1) in exp_info.iterrows():
        pm1 = float(pm1)
        node_attr = {'pepmass': pm1}
        G.add_node(id1, **node_attr)  # add nodes and attributes

    # Self clustering
    for i, (id1, pm1, charge1, spec1) in tqdm(exp_info.iterrows(), total=len(exp_info)):
        pm1 = float(pm1)
        charge1 = int(charge1)
        spec1 = spectral_entropy.clean_spectrum(spec1, max_mz=pm1 - 0.01, noise_removal = 0.01)
        mz1 = np.array(spec1[:, 0], dtype=np.float64)
        intensity1 = np.array(spec1[:, 1], dtype=np.float64)
        spectrum1 = sus.MsmsSpectrum(identifier=id1, precursor_mz=pm1 , precursor_charge= charge1
                                     , mz = mz1,intensity=intensity1).remove_precursor_peak(0.01, "Da")
        peaks1 = len(spec1)
        for j, (id2, pm2, charge2, spec2) in exp_info.iloc[:i, ].iterrows():
            pm2 = float(pm2)
            charge2 = int(charge2)
            spec2 = spectral_entropy.clean_spectrum(spec2, max_mz=pm2 - 0.01, noise_removal=0.01)
            mz2 = np.array(spec2[:, 0], dtype=np.float64)
            intensity2 = np.array(spec2[:, 1], dtype=np.float64)
            spectrum2 = sus.MsmsSpectrum(identifier=id2, precursor_mz=pm2 , precursor_charge=charge2, mz=mz2,
                                         intensity=intensity2).remove_precursor_peak(0.01, "Da")

            sim, mps, pp = 0.0, 0, 0.0
            if args.self_clustering_method == 'modified_cosine':
                try:
                    result = modified_cosine(spectrum1, spectrum2, fragment_mz_tolerance=0.02)
                    sim = result.score
                    mps = result.matches
                    pp = mps / peaks1
                except:
                    pass
            elif args.self_clustering_method == 'neutral_loss':
                try:
                    result = neutral_loss(spectrum1, spectrum2, fragment_mz_tolerance=0.02)
                    sim = result.score
                    mps = result.matches
                    pp = mps / peaks1
                except:
                    pass
            else:
                try:
                    sim = similarity(spec1, spec2, method=args.self_clustering_method, ms2_da=0.02)
                    result = modified_cosine(spectrum1, spectrum2, fragment_mz_tolerance=0.02)
                    mps = result.matches
                    pp = mps / peaks1
                except:
                    pass
            if sim >= args.self_clustering_similarity \
                    and mps >= args.self_clustering_peaks:
                edge_attr = {'pair_similarity': sim, 'matched_peaks': mps, 'peak_percentage': pp,'edge_type': args.self_clustering_method}
                G.add_edge(id1, id2, **edge_attr)
    G = mn_curating(G, args.top_k)
    print('Self clustering finished!')

    # Loading and preprocessing the in-silico MS1 match result generated by 'ms1_match()'
    isdbms1_result_path = os.path.join(parent_folder, f'IS_MS1match_{os.path.basename(args.quant_file)}')
    isms1_match_df = functions.df_preprocess(isdbms1_result_path)
    isms1_match_df['pair_similarity'] = np.nan
    isms1_match_df['mps'] = np.nan
    isms1_match_df['pp'] = np.nan

    # 0, 1, 2 represent the in-silico MS2 predicted by CFM-ID at three collision energy (10, 20 and 40 eV)
    # Pick up the
    for i in range(len(isms1_match_df)):
        max_values0 = isms1_match_df.loc[i, 'pair_similarity0']
        max_values1 = isms1_match_df.loc[i, 'pair_similarity1']
        max_values2 = isms1_match_df.loc[i, 'pair_similarity2']
        max_mps0 = isms1_match_df.loc[i, 'mps0']
        max_mps1 = isms1_match_df.loc[i, 'mps1']
        max_mps2 = isms1_match_df.loc[i, 'mps2']
        max_pp0 = isms1_match_df.loc[i, 'pp0']
        max_pp1 = isms1_match_df.loc[i, 'pp1']
        max_pp2 = isms1_match_df.loc[i, 'pp2']
        isms1_match_df.loc[i, 'pair_similarity'] = max(max_values0, max_values1, max_values2)# Filtering by maximum spectral similarity
        isms1_match_df.loc[i, 'mps'] = max(max_mps0, max_mps1, max_mps2)
        isms1_match_df.loc[i, 'pp'] = max(max_pp0, max_pp1, max_pp2)

    isms1_match_df['pair_similarity'] = pd.to_numeric(isms1_match_df['pair_similarity'], errors='coerce')# Convert missing values to N/A
    isms1_match_df['pp'] = pd.to_numeric(isms1_match_df['pp'], errors='coerce')
    index_match, index_pp_match = [], []
    index_unmatch = []

    for j in row_ids:  # traverse IS_MS1match_result by id
        temp_df = isms1_match_df[isms1_match_df['row ID'] == j]
        sim_idx = temp_df['pair_similarity'].idxmax()  # get index of maximum pair_similarity
        pp_idx = temp_df['pp'].idxmax()  # get index of maximum peak_percentage
        if not pd.isna(pp_idx):
            index_pp_match.append(pp_idx)
        if not pd.isna(sim_idx):
            index_match.append(sim_idx)
        elif pd.isna(temp_df['match_id']).all(): # get index of Feature without ms1 match
            index_unmatch.extend(temp_df.index)

    '''C2,spectral similarity annotation'''
    df_new_match = isms1_match_df.loc[index_match].reset_index(drop=True)
    df_new_match_well = df_new_match[(df_new_match['pair_similarity'] >= args.is_library_matching_similarity) & (
            df_new_match['mps'] >= args.is_library_matching_peaks)].reset_index(drop=True)
    for i in range(len(df_new_match_well)):
        pair_sim = df_new_match_well.loc[i, 'pair_similarity']
        matched_peaks = int(df_new_match_well.loc[i, 'mps'])
        peak_percentage = df_new_match_well.loc[i, 'pp']
        spec1_id = str(df_new_match_well.loc[i, 'row ID'])
        spec2_id = str(df_new_match_well.loc[i, 'match_id'])
        edge_attr = {'pair_similarity': pair_sim, 'matched_peaks': matched_peaks, 'peak_percentage': peak_percentage,
                     'edge_type': 'similarity'}
        G.add_edge(spec1_id, spec2_id, **edge_attr)
        G.nodes[spec1_id]['level'] = 'C2'
        G.nodes[spec2_id]['class'] = 'IS'
        G.nodes[spec2_id]['level'] = 'DB'
        G.nodes[spec2_id]['smile'] = df_new_match_well.loc[i, 'match_smiles']

    '''C1'''
    df_new_match_unwell = df_new_match[(df_new_match['pair_similarity'] < args.is_library_matching_similarity) | (
            df_new_match['mps'] < args.is_library_matching_peaks)].reset_index(
        drop=True)
    for i in range(len(df_new_match_unwell)):
        spec1_id = str(df_new_match_unwell.loc[i, 'row ID'])
        G.nodes[spec1_id]['level'] = 'C1'

    '''C2, peak percentage annotation'''
    df_new_pp_match = isms1_match_df.loc[index_pp_match].reset_index(drop=True)
    df_new_pp_match_well = df_new_pp_match[(df_new_pp_match['pp'] >= args.peak_percentage_threshold) & (
            df_new_pp_match['mps'] >= args.is_library_matching_peaks)].reset_index(
        drop=True)
    for i in range(len(df_new_pp_match_well)):
        pair_sim = df_new_pp_match_well.loc[i, 'pair_similarity']
        matched_peaks = int(df_new_pp_match_well.loc[i, 'mps'])
        peak_percentage = df_new_pp_match_well.loc[i, 'pp']
        spec1_id = str(df_new_pp_match_well.loc[i, 'row ID'])
        spec2_id = str(df_new_pp_match_well.loc[i, 'match_id'])
        edge_attr = {'pair_similarity': pair_sim, 'matched_peaks': matched_peaks, 'peak_percentage': peak_percentage,
                     'edge_type': 'peak_percentage'}
        G.add_edge(spec1_id, spec2_id, **edge_attr)
        G.nodes[spec1_id]['level'] = 'C2'
        G.nodes[spec2_id]['class'] = 'IS'
        G.nodes[spec2_id]['level'] = 'DB'
        G.nodes[spec2_id]['smile'] = df_new_pp_match_well.loc[i, 'match_smiles']

    # Loading and preprocessing the experimental MS1 match result generated by 'ms1_match()'
    edbms1_result_path = os.path.join(parent_folder, f'E_MS1match_{os.path.basename(args.quant_file)}')
    edbms1_match_df = functions.df_preprocess(edbms1_result_path)

    edb_index_match, edb_pp_index_match = [], []
    edb_index_unmatch = []
    edbms1_match_df['pair_similarity'] = pd.to_numeric(edbms1_match_df['pair_similarity'], errors='coerce')
    edbms1_match_df['pp'] = pd.to_numeric(edbms1_match_df['pp'], errors='coerce')
    for j in row_ids:
        temp_df = edbms1_match_df[edbms1_match_df['row ID'] == j]
        idx = temp_df['pair_similarity'].idxmax()  # get index of match with maximum pair_similarity
        pp_idx = temp_df['pp'].idxmax() # get index of match with maximum peak_percentage
        if not pd.isna(pp_idx):
            edb_pp_index_match.append(pp_idx)
        if not pd.isna(idx):
            edb_index_match.append(idx)
        elif pd.isna(temp_df['match_id']).any():  # get index of features without MS1 match
            edb_index_unmatch.extend(temp_df.index.values.tolist())

    '''B2, spectral similarity annotation'''
    edb_df_new_match = edbms1_match_df.loc[edb_index_match].reset_index(drop=True)
    edb_quant_df_new_match_well = edb_df_new_match[
        (edb_df_new_match['pair_similarity'] >= args.library_matching_similarity) & (
                edb_df_new_match['mps'] >= args.library_matching_peaks)].reset_index(drop=True)
    for i in range(len(edb_quant_df_new_match_well)):
        pair_sim = edb_quant_df_new_match_well.loc[i, 'pair_similarity']
        matched_peaks = int(edb_quant_df_new_match_well.loc[i, 'mps'])
        peak_percentage = edb_quant_df_new_match_well.loc[i, 'pp']
        spec1_id = str(edb_quant_df_new_match_well.loc[i, 'row ID'])
        spec2_id = str(edb_quant_df_new_match_well.loc[i, 'match_id'])
        edge_attr = {'pair_similarity': pair_sim, 'matched_peaks': matched_peaks, 'peak_percentage': peak_percentage,
                     'edge_type': 'similarity'}
        G.add_edge(spec1_id, spec2_id, **edge_attr)
        G.nodes[spec1_id]['level'] = 'B2'
        G.nodes[spec2_id]['class'] = 'EDB'
        G.nodes[spec2_id]['level'] = 'DB'
        G.nodes[spec2_id]['smile'] = edb_quant_df_new_match_well.loc[i, 'match_smiles']

    '''B1'''
    edb_quant_df_new_match_unwell = edb_df_new_match[
        (edb_df_new_match['pair_similarity'] < args.library_matching_similarity) | (
                edb_df_new_match['mps'] < args.library_matching_peaks)].reset_index(drop=True)
    for i in range(len(edb_quant_df_new_match_unwell)):
        spec1_id = str(edb_quant_df_new_match_unwell.loc[i, 'row ID'])
        G.nodes[spec1_id]['level'] = 'B1'

    '''B2, peak percentage annotation'''
    edb_df_new_pp_match = edbms1_match_df.loc[edb_pp_index_match].reset_index(drop=True)
    edb_quant_df_new_pp_match_well = edb_df_new_pp_match[
        (edb_df_new_pp_match['pp'] >= args.peak_percentage_threshold) & (
                edb_df_new_pp_match['mps'] >= args.library_matching_peaks)].reset_index(drop=True)
    for i in range(len(edb_quant_df_new_pp_match_well)):
        pair_sim = edb_quant_df_new_pp_match_well.loc[i, 'pair_similarity']
        matched_peaks = int(edb_quant_df_new_pp_match_well.loc[i, 'mps'])
        peak_percentage = edb_quant_df_new_pp_match_well.loc[i, 'pp']
        spec1_id = str(edb_quant_df_new_pp_match_well.loc[i, 'row ID'])
        spec2_id = str(edb_quant_df_new_pp_match_well.loc[i, 'match_id'])
        edge_attr = {'pair_similarity': pair_sim, 'matched_peaks': matched_peaks, 'peak_percentage': peak_percentage,
                     'edge_type': 'peak_percentage'}
        G.add_edge(spec1_id, spec2_id, **edge_attr)
        G.nodes[spec1_id]['level'] = 'B2'
        G.nodes[spec2_id]['class'] = 'EDB'
        G.nodes[spec2_id]['level'] = 'DB'
        G.nodes[spec2_id]['smile'] = edb_quant_df_new_pp_match_well.loc[i, 'match_smiles']

    '''Class A : MS1 no match'''
    ms1_match_file = os.path.join(parent_folder, f'MS1match_{os.path.basename(args.quant_file)}')  # MS1 result file
    ms1_match_df = functions.df_preprocess(ms1_match_file)

    temp_df = ms1_match_df.loc[:, ['row ID', 'isdbms1_id', 'edbms1_id']]
    empty_rows = temp_df[temp_df.loc[:, ['isdbms1_id', 'edbms1_id']].isnull().all(axis=1)].reset_index(drop=True)
    for i in range(len(empty_rows)):
        try:  # features in xxx_quant.csv may have no tandem mass
            spec1_id = str(empty_rows.loc[i, 'row ID'])
            G.nodes[spec1_id]['level'] = 'A'
        except:
            pass

    MN_file = os.path.join(parent_folder,
                           f'{os.path.splitext(os.path.basename(args.mgf_file))[0]}_{args.self_clustering_method}_{args.self_clustering_similarity}_{args.self_clustering_peaks}.graphml')
    nx.write_graphml(G, MN_file)
    print('Molecular networking annotation finished!')

if __name__ == '__main__':
    t = time.time()


    print(f'Finish in {(time.time() - t) / 60:.2f}min')