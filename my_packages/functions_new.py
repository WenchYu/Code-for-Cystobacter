
'''
Basic Functions for MSanalyst
'''

import os,re,ast,json,ujson
import spectral_entropy as se
import pandas as pd
import numpy as np
import spectrum_utils.spectrum as sus
from ms_entropy.file_io import spec_file
from matchms.importing import load_from_json, load_from_mgf,load_from_msp, load_from_mzml, load_from_mzxml, load_from_usi
from collections import namedtuple
from my_packages import peaktools
from tqdm import trange,tqdm
TopK = namedtuple('topk',['index','number']) # e.g. topk(index=2, number=8)

class MgfProcess:
    def __init__(self, MGF_FILE = None):
        self.MGF_FILE = MGF_FILE
        if 'GNPS' not in MGF_FILE:
            self.MGF_INFO = self.mgf_extract()
        else:
            self.MGF_INFO = self.gnps_mgf_extract()
        # if MGF_FILE is not None:
        #     self.MGF_FILE = MGF_FILE
        #     # self.MGF_FILE = self.mgf_process(MGF_FILE)
        # else:
        #     self.MGF_FILE = pd.DataFrame()

    def line_extract(self, START_TEXT):
        '''
        Extract all the lines starting with a specific keyword and save as a idlist
        :param MGF_FILE: Path including suffix of the **text** file you intend to slice
        :param START_TEXT: Starting keyword
        :return: A idlist containing content <str> after keywords
        '''
        with open(self.MGF_FILE, 'r') as f:
            CONTENT = [line[len(START_TEXT):].rstrip() for line in f if line.startswith(START_TEXT)]
        return CONTENT

    def spectra_extract(self, START_TEXT, END_TEXT, skip_words=None):
        '''
        Extracting Horizontal and vertical coordinates of tandem mass
        :param file:
        :param start_txt:
        :param end_txt:
        :param skip_words:
        :return: A idlist contain lists of sliced content, like[[],[],...,[]],and converting to an array
        '''
        if skip_words is None:
            skip_words = []
        spectra = []
        with open(self.MGF_FILE, 'r') as f:
            lines = f.readlines()
            start_idx = 0
            for i in range(len(lines)):
                if START_TEXT in lines[i]:
                    if any(word in lines[i + 1] for word in skip_words):
                        start_idx = i + 2
                    else:
                        start_idx = i + 1
                elif END_TEXT in lines[i]:
                    spectrum = ''.join(lines[start_idx:i])
                    spectra_list = spectrum.split('\n')[:-1]
                    temp = []
                    for s in spectra_list:
                        m_z, intensity = s.split()
                        temp.append([float(m_z), float(intensity)])
                    # temp = np.array(temp, dtype=np.float64) # this led to the introduction of ... and incomplete spectra
                    spectra.append(temp)
        return spectra

    def mgf_extract(self):
        '''
        Process MGF file to extract relevant information.
        :param mgf_file: '.mgf'
        :return: id<str> pepmass<str>, ms2<np array>
        '''
        ID_TEXT = 'FEATURE_ID='
        PEPMASS_TEXT = 'PEPMASS='
        CHARGE_TEXT = 'CHARGE='

        id = self.line_extract(ID_TEXT)
        pepmass = self.line_extract(PEPMASS_TEXT)
        charge = self.line_extract(CHARGE_TEXT)
        charge = [s.replace('+', '') for s in charge]

        start_txt = 'MSLEVEL=2'
        end_txt = 'END'
        ms2 = self.spectra_extract(start_txt, end_txt, skip_words=['MERGED'])

        exp_info = pd.DataFrame({
            'id': id,
            'pepmass': pepmass,
            'charge': charge,
            'ms2': ms2
        })
        exp_info = exp_info[exp_info['ms2'].apply(len) > 1]  # delete empty idlist
        exp_info = exp_info.reset_index(drop=True)  # reindex
        return exp_info

    def feature_extract(self, FEATURE_ID):
        '''
        Retrieve information from MGF file based on ID.
        :param mgf_id:
        :return: pepmass<float>, spec<np.array>, spectrum<vsl object>
        '''
        mgf_info = self.MGF_INFO
        try:
            pepmass = float(mgf_info[mgf_info['id'] == FEATURE_ID]['pepmass'].iloc[0])
            charge = int(mgf_info[mgf_info['id'] == FEATURE_ID]['charge'].iloc[0])
            spec = mgf_info[mgf_info['id'] == FEATURE_ID]['ms2'].iloc[0]
            spec = se.clean_spectrum(spec, max_mz=pepmass + 0.01)

            mz = np.array(spec[:, 0])
            intensity = np.array(spec[:, 1])
            spectrum = sus.MsmsSpectrum(
                identifier=FEATURE_ID,
                precursor_mz=pepmass,
                precursor_charge=charge,
                mz=mz,
                intensity=intensity
            )
            return {
                'pepmass': pepmass,
                'spec': spec,
                'spectrum': spectrum,
                'charge': charge,
                'id': FEATURE_ID
            }
        except:
            raise ValueError(f"No data found for mgf_id: {FEATURE_ID}")

    def gnps_mgf_extract(self):
        '''
        Process MGF file to extract relevant information.
        :param mgf_file: '.mgf'
        :return: a dict used to generate json file
        id<str> pepmass<str>, ms2<np array>
        '''

        SPECTRUMID = self.line_extract('SPECTRUMID=')
        PEPMASS = self.line_extract('PEPMASS=')
        CHARGE = self.line_extract('CHARGE=')
        MSLEVEL = self.line_extract('MSLEVEL=')
        IONMODE = self.line_extract('IONMODE=')
        NAME = self.line_extract('NAME=')
        SMILE = self.line_extract('SMILES=')
        INSTRUMENT = self.line_extract('SOURCE_INSTRUMENT=')
        # charge = [s.replace('+', '') for s in charge]

        START_TEXT = 'SCAN'
        END_TEXT = 'END'
        MS2_SPECTRUM = self.spectra_extract(START_TEXT, END_TEXT, skip_words=['MERGED'])

        dict = {}
        print('Start to convert to dict')
        for i in trange(len(SPECTRUMID)):
            dict[f'{SPECTRUMID[i]}'] = {
                'PEPMASS': f'{PEPMASS[i]}',
                'CHARGE': f'{CHARGE[i]}',
                'MSLEVEL': f'{MSLEVEL[i]}',
                'IONMODE': f'{IONMODE[i]}',
                'NAME': f'{NAME[i]}',
                'SMILE': f'{SMILE[i]}',
                'INSTRUMENT': f'{INSTRUMENT[i]}',
                'MS2_SPECTRUM': f'{MS2_SPECTRUM[i]}'
                }

        return dict

    def query_mass_process(self, QMS1, QMS2):
        '''
        Specifically designed for MSanalyst
        Directly process input query MS1 and MS2 spectra
        :param qms1: e.g. '381.2958'
        :param qms2: e.g. '381.2284 1.0E2 381.2344 1.1E2 381.2822 1.1E2 381.2842 1.3E2 381.2862 5.2E2'
        :return: e.g. '381.2284 1.0E2 381.2344 1.1E2 381.2822 1.1E2 381.2842 1.3E2 381.2862 5.2E2'
        '''
        id = '1'
        pepmass = QMS1
        charge = '1'
        try:
            spectra = []
            temp = []
            lines = QMS2.strip().split('\n')
            for line in lines:
                m_z, intensity = line.split()
                temp.append([float(m_z), float(intensity)])
            temp = np.array(temp, dtype=np.float64)
            spectra.append(temp)
        except:
            spectra = []
            temp = []
            elements = QMS2.split()
            for i in range(0, len(elements), 2):
                temp.append([float(elements[i]), float(elements[i + 1])])
            temp = np.array(temp, dtype=np.float64)
            spectra.append(temp)

        exp_info = pd.DataFrame({
            'id': [id],
            'pepmass': [pepmass],
            'charge': [charge],
            'ms2': [spectra[0]]  # spectra is a idlist containing one numpy array
        })
        return exp_info

def load_spectra_from_file(INPUT_FILE):
    '''
    It will return a list containing spectral objects
    Properties can be accessed by (.mz, .intensities, .metadata)
    '''
    try:
        FILE_TYPE = spec_file.guess_file_type_from_file_name(INPUT_FILE) # msp, mgf, mzml, mzml, hdf5, raw, lbm2
    except:
        FILE_TYPE = INPUT_FILE
    if FILE_TYPE == 'json':
        SPECTRA = list(load_from_json(INPUT_FILE))
    if FILE_TYPE == 'mgf':
        SPECTRA = list(load_from_mgf(INPUT_FILE))
    if FILE_TYPE == 'msp':
        SPECTRA = list(load_from_msp(INPUT_FILE))
    if FILE_TYPE == 'mzml':
        SPECTRA = list(load_from_mzml(INPUT_FILE))
    if FILE_TYPE == 'mzxml':
        SPECTRA = list(load_from_mzxml(INPUT_FILE))
    if 'mzspec:' in INPUT_FILE:
         SPECTRA = list(load_from_usi(INPUT_FILE))
    return SPECTRA

def spec_str2array(SPEC_STR,PEPMASS):
    '''
    Convert MS2_spectrum_str to MS2_spectrum_array
    :param SPEC_STR:
    :param PEPMASS:
    :return:
    '''
    SPEC_STR = SPEC_STR.replace(',', ' ').replace('[', ' ').replace(']', ' ')  # Formatting the spec_str[mz ins mz ins mz ins ...]
    SPEC_STR = SPEC_STR.split()  # Split by space
    SPEC_FLOAT = [float(x) for x in SPEC_STR]  # floating
    SPEC_ARRAY = np.array(SPEC_FLOAT).reshape(-1, 2)  # mz ins\n mz ins\n ... <np.array>
    SPEC_ARRAY = se.clean_spectrum(SPEC_ARRAY, max_mz=PEPMASS + 0.01)
    return SPEC_ARRAY

def GNPS_info_format(GNPS_INFO,CCMSID):
    '''
    Json only support text format, conversion to np.array is needed.
    :param GNPS_INFO:
    :param CCMSID:
    :return:
    '''
    EX_INFO = GNPS_INFO[CCMSID]
    spec_str = EX_INFO['MS2_SPECTRUM']
    PEPMASS = float(EX_INFO['PEPMASS'])

    spec_str = spec_str.replace(',', '').replace('[', '').replace(']', '')  # Formatting the spec_str[mz ins mz ins mz ins ...]
    spec_str = spec_str.split()  # Split by space
    spec_float = [float(x) for x in spec_str]  # floating
    spec_array = np.array(spec_float).reshape(-1, 2)  # mz ins\n mz ins\n ... <np.array>
    spec_array = se.clean_spectrum(spec_array, max_mz=PEPMASS + 0.01)
    mz = np.array(spec_array[:, 0])
    intensity = np.array(spec_array[:, 1])

    SPECTRUM= sus.MsmsSpectrum(
        identifier=CCMSID,
        precursor_mz=PEPMASS,
        precursor_charge=1,
        mz=mz,
        intensity=intensity
    )
    return spec_array,SPECTRUM

def json_load(JSONFILE):
    with open(JSONFILE,'r') as f:
        return ujson.load(f)

def refms_compare(GNPS_INFO, ALGORITHM, CCMSID1, CCMSID2):

    SPEC1, SPECTRUM1 = GNPS_info_format(GNPS_INFO, CCMSID1)
    SPEC2, SPECTRUM2 = GNPS_info_format(GNPS_INFO, CCMSID2)

    if ALGORITHM == 'cosine':
        return peaktools.cosine(SPECTRUM1, SPECTRUM2, 0.05)
    elif ALGORITHM == 'modified_cosine':
        return peaktools.modified_cosine(SPECTRUM1, SPECTRUM2, 0.05)
    elif ALGORITHM == 'neutral_loss':
        return peaktools.neutral_loss(SPECTRUM1, SPECTRUM2, 0.05)
    else:
        return se.similarity(SPEC1, SPEC2, method=ALGORITHM, ms2_da=0.05)


def arrary2list(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {obj.__class__.__name__} is not JSON serializable")

def create_result_folders(args):
    '''
    Create result folders based on the quant file
    '''
    df = pd.read_csv(args.quant_file)
    parent_folder = f'{args.output}/{os.path.splitext(os.path.basename(args.quant_file))[0]}_result'# output/_quant_result/**
    os.makedirs(parent_folder, exist_ok=True)
    for _, row in df.iterrows():
        folder_name = f"{parent_folder}/{int(row['row ID'])}"
        os.makedirs(folder_name, exist_ok=True)
    print('Result folders have been created!')

def create_subresults(args):
    '''
    Split the results  after the MS1 match
    Create a separate CSV for each row ID, writing the corresponding information to facilitate detailed inspection
    '''
    parent_folder =  f'{args.output}/{os.path.splitext(os.path.basename(args.quant_file))[0]}_result' # filename ''output/_quant_result/**''
    npms1_result_path =os.path.join(parent_folder, f'IS_MS1match_{os.path.basename(args.quant_file)}')
    edbms1_result_path = os.path.join(parent_folder, f'E_MS1match_{os.path.basename(args.quant_file)}')

    quant_df = df_preprocess(args.quant_file)
    npms1_match_df = df_preprocess(npms1_result_path)
    edbms1_match_df = df_preprocess(edbms1_result_path)


    for i in range(len(quant_df)):
        id = quant_df['row ID'][i]
        folder_name = os.path.join(parent_folder, str(id))

        npcsv_file = os.path.join(folder_name, f'IS_MS1match_{str(id)}.csv') # isdb results
        if not os.path.exists(npcsv_file):
            pd.DataFrame(columns=npms1_match_df.columns).to_csv(npcsv_file, index=False)
        selected_rows =npms1_match_df.loc[npms1_match_df['row ID'] == id]
        with open(npcsv_file, 'a', newline='') as f1:
            selected_rows.to_csv(f1, index=False, header=False)

        edbcsv_file = os.path.join(folder_name, f'E_MS1match_{str(id)}.csv') # edb result
        if not os.path.exists(edbcsv_file):
            pd.DataFrame(columns=edbms1_match_df.columns).to_csv(edbcsv_file, index=False)
        selected_rows = edbms1_match_df.loc[edbms1_match_df['row ID'] == id]
        with open(edbcsv_file, 'a', newline='') as f2:
            selected_rows.to_csv(f2, index=False, header=False)

def get_edb_info(gnps_info, gnps_id):
    '''

    :param isdb_info:
    :param id:
    :return:
    '''
    keys_to_retrieve = ['smiles', 'pepmass', 'ms2','charge']
    values = [gnps_info[gnps_id][key] for key in keys_to_retrieve]
    smiles, pepmass, spec, charge = values
    # string convertion
    pepmass = float(pepmass)
    charge = int(charge)
    spec = np.asarray(ast.literal_eval(spec))
    mz = np.array(spec[:, 0])
    spectrum = sus.MsmsSpectrum(identifier=f'{gnps_id}'
                                 , precursor_mz=pepmass
                                 , precursor_charge=charge
                                 , mz=mz
                                 , intensity=spec[:, 1])

    return {'smiles': smiles, 'pepmass': pepmass
        , 'spec': spec, 'spectrum': spectrum,'charge': charge}

def get_isdb_info(isdb_info, is_id):
    '''

    :param isdb_info:
    :param id:
    :return:
    '''
    keys_to_retrieve = ['smiles', 'pepmass', 'energy0_ms2', 'energy1_ms2', 'energy2_ms2']
    values = [isdb_info[is_id][key] for key in keys_to_retrieve]
    smiles, pepmass, e0spec, e1spec, e2spec = values
    # string convertion
    pepmass = float(pepmass)
    e0spec = np.asarray(ast.literal_eval(e0spec))
    e1spec = np.asarray(ast.literal_eval(e1spec))
    e2spec = np.asarray(ast.literal_eval(e2spec))

    mz0 = np.array(e0spec[:, 0])
    spectrum0 = sus.MsmsSpectrum(identifier=f'e0_{is_id}'
                                 , precursor_mz=pepmass
                                 , precursor_charge=1
                                 , mz=mz0
                                 , intensity=e0spec[:, 1])
    mz1 = np.array(e1spec[:, 0])
    spectrum1 = sus.MsmsSpectrum(identifier = f'e1_{is_id}'
                                 , precursor_mz=pepmass
                                 , precursor_charge=1
                                 , mz=mz1
                                 , intensity=e1spec[:, 1])
    mz2 = np.array(e2spec[:, 0])
    spectrum2 = sus.MsmsSpectrum(identifier=f'e2_{is_id}'
                                 , precursor_mz=pepmass
                                 , precursor_charge=1
                                 , mz=mz2
                                 , intensity=e2spec[:, 1])

    return {'smiles': smiles, 'pepmass': pepmass
        , 'e0spec': e0spec, 'e1spec': e1spec, 'e2spec': e2spec
        , 'e0spectrum': spectrum0, 'e1spectrum': spectrum1, 'e2spectrum': spectrum2}

def df_preprocess(filename):
    '''
    Preprocess DataFrame by removing empty columns and resetting index.
    '''
    if filename.endswith('.csv'):
        df = pd.read_csv(filename, low_memory=False)
    elif filename.endswith('.tsv'):
        df = pd.read_csv(filename, sep='\t', low_memory=False)
    elif filename.endswith('.xlsx') or filename.endswith('.xls'):
        df = pd.read_excel(filename)
    else:
        raise ValueError("Unsupported file format. Please use .csv, .tsv, or .xlsx files.")

    if  df.index[-1] != len(df)-1:
        df.index.name = ''
        df.reset_index(inplace=True)
    return df

def calculate_ppm(query_mass_value: float, reference_mass_value: float) -> float:
    '''
    Calculate parts per million (ppm) for mass values.
    '''
    if not isinstance(query_mass_value, (int, float)) or not isinstance(reference_mass_value, (int, float)):
        raise TypeError('Input parameters must be numbers.')
    if reference_mass_value != 0:
        return abs((query_mass_value - reference_mass_value) / reference_mass_value * 1e6)
    return float('inf')

def db_parsing():
    '''
    Parse default databases of MSanalyst.
    '''
    isdb_file = './msdb/isdb_info.json'
    edb_file = './msdb/edb_info.json'
    with open(isdb_file, 'r') as f:
        isdb_info = json.load(f)
    with open(edb_file, 'r') as f1:
        gnps_info = json.load(f1)
    return isdb_info, gnps_info

def list_files(directory,keyword):
    '''
    idlist files with keyword
    :param directory: dirt
    :return: A idlist containing all .graphml file paths
    '''
    graphml_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(keyword):
                graphml_files.append(os.path.join(root, file))
    return graphml_files

def ex_algorithm_name(file,prefix):
    '''
    extract algorithm name from '.graphml' files
    e.g. dot_product from KutzOsmac_dot_product_0.7_3.graphml
    :return:
    '''
    match = re.search(r'^(\w+)_|(\w+)_\d+\.\d+', file)
    if match:
        pattern = match.group(1) if match.group(1) else match.group(2)
        pattern = pattern.replace(f'{prefix}_', '')
        return pattern


if __name__ == '__main__':
    test_mgf = '../notebooks/data/kutz/KutzOsmac.mgf'
    mgf_processor = MgfProcess(test_mgf)

    print('')