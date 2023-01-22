import pandas as pd
import numpy as np

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

def get_pos_ctrl_corrected(corrected_sample, corrected_pos_ctrl):
    
    return corrected_sample.div(corrected_pos_ctrl.mean(axis=1), axis=0)

def read_map(plate_map):

    '''serialize a map file into a dataframe'''
    plate_map = plate_map.melt(id_vars=['group']).dropna()
    plate_map['well'] = plate_map['group'] + plate_map['variable'].astype(int).apply(lambda x: '{:02d}'.format(x))
    return dict(zip(plate_map['well'], plate_map['value']))

def setup_header(raw_data):

    '''check if the data has single or multiple group'''
    start_idx = 3 if 'Group' in raw_data.columns else 2

    '''check if the data has single or double header (single vs time-series measurement)'''
    if pd.isna(raw_data.iloc[0, 0]) or raw_data.iloc[0, 0]!='A01':
        raw_data.columns = raw_data.columns.tolist()[:start_idx] + raw_data.iloc[0].tolist()[start_idx:]
        raw_data.drop(0, inplace=True)
    
    return raw_data

def map_column_name(data, sample_map):
    
    data['name'] = data['Well'].apply(lambda x: sample_map[x] if x in sample_map else x)
    data = data.set_index('name').drop('Well', axis=1)
    return data

def parse_data(data, num_data=4):
    
    data = data.drop(['Content'], axis=1)
    
    datas = []
    idxs = [0]
    timepoint = int(data.shape[1]/num_data)
    
    for i in range(1, num_data+1):
        idxs.append(idxs[i-1] + timepoint)
        datas.append((data.iloc[:, idxs[i-1]:idxs[i]]).astype(float).T)

    return datas

def parse_per_content(data):

    all_data = {}
    
    for content in list(set(data['Content'].apply(lambda x: x.split()[0]))):
        all_data[content] = parse_data(data[data['Content'].str.startswith(content)])

    return all_data

def get_chunked_data(data):
    
    chunked_data = {}
    
    if 'Group' not in data.columns:
        data['Group'] = 'A'
        
    for group in data['Group'].unique():
        chunked_data[group] = parse_per_content(data[(data['Group']==group)].drop('Group', axis=1))
        
    return chunked_data

def generate_raw(folder, filename='raw', mapname='plate_map'):
    
    sample_map = read_map(pd.read_csv('datasets/experiment/{}/{}.csv'.format(folder, mapname)))
    raw_data = setup_header(pd.read_csv('datasets/experiment/{}/{}.csv'.format(folder, filename)))
    mapped_data = map_column_name(raw_data, sample_map)
    chunked_data = get_chunked_data(mapped_data)

    groups = list(chunked_data.keys())
    for g in groups:
        keys = list(chunked_data[g].keys())
        print(g, keys)
        
    return chunked_data

def generate_neg_corrected(raw_data, sample_group=['A'], pos_group='B', fluo_idx=3, od_idx=0):

    samples = []
    for s in sample_group:
        data = raw_data[s]
        samples.append(get_neg_ctrl_corrected([data['Sample'][fluo_idx], data['Blank'][fluo_idx], data['Negative'][fluo_idx],
                data['Sample'][od_idx], data['Blank'][od_idx], data['Negative'][od_idx]]))
    
    if pos_group is not None:
        data = raw_data[pos_group]
        pos_control = get_neg_ctrl_corrected([data['Positive'][fluo_idx], data['Blank'][fluo_idx], data['Negative'][fluo_idx],
                data['Positive'][od_idx], data['Blank'][od_idx], data['Negative'][od_idx]])
    else:
        pos_control = None

    return samples, pos_control

def generate_pos_corrected(samples, pos_control):

    if pos_control is None:
        print('WARNING: there is no positive control!')
        return samples

    corrected_data = []
    for sample in samples:
        corrected_data.append(get_pos_ctrl_corrected(sample, pos_control))
    return corrected_data

def get_data_at(data, h=8, per_minute=20):
    
    timepoint = int((60/per_minute) * h)
    return np.abs(data.iloc[timepoint])


def inverse_hill(rpu, ag, K, n, eps):
    
    ag_, K_, n_, eps_ = 10**ag, 10**K, 10**n, 10**eps
    #return (((rpu * (K_**n_)) - (ag_ * eps_ * (K_**n_))) / (ag_ - rpu))**(1/n_)
    return (((rpu - (ag_ * eps_)) * (K_**n_))/(ag_ - rpu))**(1/n_)