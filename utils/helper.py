import pandas as pd
import numpy as np
import string

### BLANK-CORRECTED AND NEGATIVE-CONTROL-CORRECTED DATA ###

def get_blank_corrected(data, blank_id='blank'):
    
    blank_corrected = data.sub(data[filter(lambda x: x.startswith(blank_id), data.columns)].mean(axis=1), axis=0)
    return blank_corrected.drop(filter(lambda x: x.startswith(blank_id), blank_corrected.columns), axis=1).reset_index(drop=True)

def get_neg_ctrl_corrected(fluo, od, neg_ctrl_id='negative-control'):
    
    blank_corr_fluo = get_blank_corrected(fluo)
    blank_corr_od = get_blank_corrected(od)
    data = blank_corr_fluo / blank_corr_od
    neg_corrected = data.sub(data[filter(lambda x: x.startswith(neg_ctrl_id), data.columns)].mean(axis=1), axis=0)
    return neg_corrected.drop(filter(lambda x: x.startswith(neg_ctrl_id), neg_corrected.columns), axis=1)

### DATA PREPROCESSING ###

def setup_header(raw_data, start_idx=3):

    '''check if the data has a single header (i.e. single measurement)'''
    #if the first row first col starts with a well, then it is a single measurement data
    if len(list(filter(raw_data.iloc[0, 0].startswith, list(string.ascii_uppercase)[:8])))==0:
        raw_data.columns = raw_data.columns.tolist()[:start_idx] + raw_data.iloc[0].tolist()[start_idx:]
        raw_data.drop(0, inplace=True)
    return raw_data

def read_map(plate_map):
    '''serialize a map file into a dataframe'''
    
    plate_map = plate_map.melt(id_vars=['group'])
    plate_map['variable'] = plate_map['variable'].astype(int)
    plate_map['Well'] = plate_map['group'] + plate_map['variable'].apply(lambda x: '{:02d}'.format(x))
    return plate_map[['Well', 'value']].dropna()

def read_dict(dictionary):
    
    sample_map = pd.Series(dictionary['short_name'].values, index=dictionary['id']).to_dict()
    sample_map.update({
        'BK': 'blank-kan',
        'BA': 'blank-amp',
        '3K3-N': 'negative-control-kan',
        '4A3-N': 'negative-control-amp',
        '4A3-P': 'positive-control-amp'
    })
    return sample_map

def generate_metadata(well, plate_map, sample_map):
    
    df = pd.merge(well, plate_map, on='Well', how='left').dropna(subset=['value']).reset_index(drop=True)
    df['short_name'] = df['value'].map(sample_map) + '_' + df['suffix'].astype(str)
    return df.dropna()

def transpose_data(df, col):
    
    df.set_index(col, inplace=True)
    df = df.transpose().reset_index()
    df = df.set_index('index')
    return df

def generate_data(df, name, datapoint=2, num_data=4, start_idx=3, col='short_name'):
    
    datas = []
    idxs = [start_idx]
    for i in range(1, num_data+1):
        idxs.append(idxs[i-1] + datapoint)
        data = (df.iloc[:, idxs[i-1]:idxs[i]]).astype(float)
        data = pd.concat([name, data], axis=1)
        data = transpose_data(data, col)
        datas.append(data)
    return datas