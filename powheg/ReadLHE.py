'''
Script with utils methods to manage les houche files
'''

import pandas as pd
import numpy as np

def ReadLHEFile(fileName='', nRowsToSkip=49):
    '''
    Helper method to read lhe files with pandas
    
    Parameters
    -------------
    - fileName: name of les houche event file wuth full path
    - nRowsToSkip: number of initial rows to skip
    
    Returns
    -------------
    - df: dataframe with events 
    '''

    print(f'\nRead input file {fileName}')

    df = pd.read_csv(fileName, skiprows=nRowsToSkip, comment='<', delim_whitespace=True, error_bad_lines=False, low_memory=False, header=0,
                     names=['pdg', 'unused0', 'unused1', 'unused2', 'unused3', 'unused4', 'px', 'py', 'pz', 'E', 'm', 'unused5', 'spin'])
    df = df[~df['pdg'].str.contains('#')]

    print('Attaching event number...', end='')
    evNum = []
    iEv = 0
    for spin in df['spin'].values:
        if np.isnan(spin):
            iEv += 1
        evNum.append(iEv)
    print('\rAttaching event number [Done]\n\n', end='', flush=True)

    df['event'] = pd.Series(evNum)
    df = df.dropna()
    df = df[['pdg', 'px', 'py', 'pz', 'E', 'm', 'spin', 'event']]

    df = df.astype({'pdg': 'int32'})

    return df
