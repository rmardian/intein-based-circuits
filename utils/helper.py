import pandas as pd
import numpy as np
import string

### BLANK-CORRECTED AND NEGATIVE-CONTROL-CORRECTED DATA ###

def get_blank_corrected(data, blank):
    
    return data.sub(blank.mean(axis=1), axis=0)
    
def get_neg_ctrl_corrected(sample):

    sample_fluo, blank_fluo, neg_fluo, sample_od, blank_od, neg_od = sample
    blank_corr_sample_fluo = get_blank_corrected(sample_fluo, blank_fluo)
    blank_corr_sample_od = get_blank_corrected(sample_od, blank_od)
    blank_corr_neg_fluo = get_blank_corrected(neg_fluo, blank_fluo)
    blank_corr_neg_od = get_blank_corrected(neg_od, blank_od)
    blank_corr_sample = blank_corr_sample_fluo / blank_corr_sample_od
    blank_corr_neg = blank_corr_neg_fluo / blank_corr_neg_od
    return blank_corr_sample.sub(blank_corr_neg.mean(axis=1), axis=0)

def get_pos_ctrl_corrected(sample, pos_ctrl):
    
    corrected_sample = get_neg_ctrl_corrected(sample)
    corrected_pos_ctrl = get_neg_ctrl_corrected(pos_ctrl)
    return corrected_sample.div(corrected_pos_ctrl.mean(axis=1), axis=0)

### DATA PREPROCESSING ###

def setup_header(raw_data):

    '''check if the data has a single or multiple group'''
    start_idx = 3 if 'Group' in raw_data.columns else 2

    '''check if the data has a single or double header (single vs time-series measurement)'''
    if pd.isna(raw_data.iloc[0, 0]):
        raw_data.columns = raw_data.columns.tolist()[:start_idx] + raw_data.iloc[0].tolist()[start_idx:]
        raw_data.drop(0, inplace=True)
    
    return raw_data

def read_map(plate_map, dictionary_file=None):

    '''serialize a map file into a dataframe'''
    plate_map = plate_map.melt(id_vars=['group'])
    plate_map['variable'] = plate_map['variable'].astype(int)
    plate_map['Well'] = plate_map['group'] + plate_map['variable'].apply(lambda x: '{:02d}'.format(x))
    plate_map.dropna(inplace=True)

    if dictionary_file is not None:
        dictionary = dict(zip(dictionary_file['id'], dictionary_file['short_name']))
        plate_map['value'] = plate_map['value'].apply(lambda x: dictionary[x] if x in dictionary else x)

    return dict(zip(plate_map['Well'], plate_map['value']))

def parse_data(data, sample_map, content='Sample', kind='control', num_data=4):

    if len(data)==0:
        return []

    data['sample'] = data['Well'].map(sample_map)
    data = data.dropna().reset_index(drop=True)
    data['row'] = data['Well'].str[0]
    data['col'] = data['Well'].str[1:].astype(int)

    if content=='Sample' and kind=='3-input on-off':
        data['suffix'] = (data['row'].apply(lambda x: str(ord(x)-65)))
    elif content=='Sample' and kind=='2-input on-off':
        data['suffix'] = (data['row'].apply(lambda x: str(ord(x)-65))%4)
    elif content=='Sample' and kind=='3-input induction matrix':
        '''4-induction level for all 3 inputs'''
        induction_states = []
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    induction_states.append(str(i) + str(j) + str(k))
        data['suffix'] = induction_states
    elif content=='Sample' and kind=='2-input induction matrix':
        '''6-induction level for both inputs'''
        data['suffix'] = (data['row'].apply(lambda x: str(ord(x)-65))) + ((data['col']-1)%6).astype(str)
    elif content=='Sample' and kind=='1-input response function':
        data['suffix'] = ((data['col']-1)%12).astype(str)
    else:
        data['suffix'] = np.arange(len(data)).astype(str)

    data['name'] = data['sample'] + '_' + data['suffix']
    data.set_index('name', inplace=True)
    if 'Group' in data.columns:
        data.drop(['Well', 'Content', 'Group', 'sample', 'row', 'col', 'suffix'], axis=1, inplace=True)
    else:
        data.drop(['Well', 'Content', 'sample', 'row', 'col', 'suffix'], axis=1, inplace=True)

    #transposing data
    datas = []
    idxs = [0]
    timepoint = int(data.shape[1]/num_data)
    for i in range(1, num_data+1):
        idxs.append(idxs[i-1] + timepoint)
        datas.append((data.iloc[:, idxs[i-1]:idxs[i]]).astype(float).T)

    return datas

def generate_data(raw_data, sample_map, kind, fluo_idx=3, od_idx=0):

    all_data = []
    chunked_data = []

    if 'Group' in raw_data.columns:
        for group in raw_data['Group'].unique():
            chunked_data.append((group, raw_data[(raw_data['Group']==group)]))
    else:
        chunked_data = [('A', raw_data.copy())]

    for group, data in chunked_data:
        
        samples = parse_data(data[data['Content'].str.startswith('Sample')],
                            sample_map, 'Sample', kind=kind)
        blanks = parse_data(data[data['Content'].str.startswith('Blank')],
                            sample_map, 'Blank')
        negatives = parse_data(data[data['Content'].str.startswith('Negative')],
                            sample_map, 'Negative')
        positives = parse_data(data[data['Content'].str.startswith('Positive')],
                            sample_map, 'Positive')

        if len(samples)>0 and len(blanks)>0 and len(negatives)>0:

            all_data.append((group, 'sample', get_neg_ctrl_corrected([samples[fluo_idx], blanks[fluo_idx], negatives[fluo_idx],
                                    samples[od_idx], blanks[od_idx], negatives[od_idx]])))

        if len(positives)>0 and len(blanks)>0 and len(negatives)>0:

            all_data.append((group, 'positive', get_neg_ctrl_corrected([positives[fluo_idx], blanks[fluo_idx], negatives[fluo_idx],
                                    positives[od_idx], blanks[od_idx], negatives[od_idx]])))

    return all_data

def generate_data_at_t_2(folder, filename, sample_map_, kind='1-input response function', pos_ctrl_map=[(0, 1)], h=8, rpu=True, single_timepoint=False):
    
    timepoint = 3*h

    raw_data = setup_header(pd.read_csv('datasets/experiment/{}/{}.csv'.format(folder, filename)))
    sample_map = read_map(sample_map_,
                          pd.read_csv('datasets/dictionary.csv'))
    data = generate_data(raw_data, sample_map, kind)

    #print([(a[0], a[1]) for a in data])
    
    ctrl_corrected = []
    for c in pos_ctrl_map:
        if rpu:
            samples = data[c[0]][2].div(data[c[1]][2].mean(axis=1), axis=0)
        else:
            samples = data[c[0]][2].copy()
        
        if not single_timepoint:
            ctrl_corrected.append(np.abs(samples.iloc[timepoint]))
        else:
            ctrl_corrected.append(np.abs(samples.mean()))
        
    return ctrl_corrected



def generate_data_at_t(folder, filename, map_file, kind='1-input response function', pos_ctrl_map=[(0, 1)], h=8, rpu=True, single_timepoint=False):
    
    timepoint = 3*h

    raw_data = setup_header(pd.read_csv('datasets/experiment/{}/{}.csv'.format(folder, filename)))
    sample_map = read_map(pd.read_csv('datasets/experiment/{}/{}.csv'.format(folder, map_file)),
                          pd.read_csv('datasets/dictionary.csv'))
    data = generate_data(raw_data, sample_map, kind)

    print([(a[0], a[1]) for a in data])
    
    ctrl_corrected = []
    for c in pos_ctrl_map:
        if rpu:
            samples = data[c[0]][2].div(data[c[1]][2].mean(axis=1), axis=0)
        else:
            samples = data[c[0]][2].copy()
        
        if not single_timepoint:
            ctrl_corrected.append(np.abs(samples.iloc[timepoint]))
        else:
            ctrl_corrected.append(np.abs(samples.mean()))
        
    return ctrl_corrected

### METRICS ###

def sum_squared_error(y_test, y_pred):

    squared_errors = [(yt - yp) ** 2 for yt, yp in zip(y_test, y_pred)]
    return np.sum(squared_errors)