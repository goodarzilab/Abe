import numpy as np
import pandas as pd
# from time import time
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, spearmanr
import plotly as py
import plotly.tools as tls
from cohort import make_study_dataset


def get_SLdataset(SL_thr = 3):
    """
    Read data from Luke's paper "Mapping the Genetic Landscape of Human Cells"
    https://doi.org/10.1016/j.cell.2018.06.010
    Table S5
    """
    myfile = ("CRISPRi_Mapping_paper/Table_S5.xlsx")
    xl = pd.ExcelFile(myfile)
    sheets = {sheet: xl.parse(sheet) for sheet in xl.sheet_names}
    # print ('Sheets:')
    # print (xl.sheet_names)
    # Read data from *gene GI scores sheet*:
    GIsheet = sheets['gene GI scores and correlations']
    raw_data = pd.DataFrame(data = {
        'Gene1': np.array(GIsheet['Unnamed: 0'][3:]),
        'Gene2': np.array(GIsheet['Unnamed: 1'][3:]),
        'K562': np.array(GIsheet['K562.4'][3:]),                # g1 <-> g2 GI scores
        'Jurkat': np.array(GIsheet['Jurkat.4'][3:])             # g1 <-> g2 GI scores
    })
    raw_K562   = raw_data.drop(columns='Jurkat').dropna()
    raw_Jurkat = raw_data.drop(columns='K562').dropna()
    K562   = raw_K562[(raw_K562['K562'] > SL_thr) &
                      (raw_K562['Gene1'] != raw_K562['Gene2'])].reset_index(drop=True)
    Jurkat = raw_Jurkat[(raw_Jurkat['Jurkat'] > SL_thr) &
                        (raw_Jurkat['Gene1'] != raw_Jurkat['Gene2'])].reset_index(drop=True)
    print (f'K562: {round(100* len(K562) / len(raw_K562), 2)}%', end = '\t')
    print (f'{len(K562)} SLs from {len(raw_K562)} unique gene pairs')
    print (f'Jurkat: {round(100* len(Jurkat) / len(raw_Jurkat), 2)}%', end = '\t')
    print (f'{len(Jurkat)} SLs from {len(raw_Jurkat)} unique gene pairs ')
    return {'K562':K562, 'Jurkat':Jurkat}


def SL_stat_test(filename, cell_line='K562'):
    SLdata = get_SLdataset()
    RNASeq, data = make_study_dataset(filename, get_RNASeq=True)
    tmp = SLdata[cell_line]
    # filter out missing gene names between GI and expression studies
    tmp = tmp[pd.DataFrame(tmp.Gene1.tolist()).isin(data.index.tolist()).any(1)].reset_index(drop=True)
    tmp = tmp[pd.DataFrame(tmp.Gene2.tolist()).isin(data.index.tolist()).any(1)].reset_index(drop=True)
    G1 = tmp.Gene1;
    G2 = tmp.Gene2;
    data_G1 = np.array(data.loc[G1], dtype=int)
    data_G2 = np.array(data.loc[G2], dtype=int)
    Cors = [spearmanr(g1, g2)[0] for g1, g2 in zip(np.array(RNASeq.loc[G1]), np.array(RNASeq.loc[G2]))]
    test_res = pd.DataFrame(data=np.stack((
        G1, G2, tmp[cell_line], Cors,
        np.sum((data_G1 < 0) * (data_G2 < 0), axis=1),  # low_low
        np.sum((data_G1 > 0) * (data_G2 < 0), axis=1),  # high_low
        np.sum((data_G1 < 0) * (data_G2 > 0), axis=1),  # low_high
        np.sum((data_G1 > 0) * (data_G2 > 0), axis=1)),  # high_high
        axis=1), columns=[
        'Gene1', 'Gene2', 'GI_Score', 'Correlation',
        'Obs_Low_Low', 'Obs_High_Low', 'Obs_Low_High', 'Obs_High_High'
        #         , 'ChiPval', 'Chi2'
    ])

    # random gene pairs
    r1, r2 = np.random.randint(0, high=len(data) - 1, size=(2, len(tmp)), dtype='int')
    rG1 = [data.index.tolist()[i] for i in r1]
    rG2 = [data.index.tolist()[i] for i in r2]
    data_rG1 = np.array(data.loc[rG1], dtype=int)
    data_rG2 = np.array(data.loc[rG2], dtype=int)
    rCors = [spearmanr(g1, g2)[0] for g1, g2 in zip(np.array(RNASeq.loc[rG1]), np.array(RNASeq.loc[rG2]))]
    random_test_res = pd.DataFrame(data=np.stack((
        rG1, rG2, rCors,
        np.sum((data_rG1 < 0) * (data_rG2 < 0), axis=1),  # low_low
        np.sum((data_rG1 > 0) * (data_rG2 < 0), axis=1),  # high_low
        np.sum((data_rG1 < 0) * (data_rG2 > 0), axis=1),  # low_high
        np.sum((data_rG1 > 0) * (data_rG2 > 0), axis=1)),  # high_high
        axis=1), columns=[
        'Gene1', 'Gene2', 'Correlation',
        'Obs_Low_Low', 'Obs_High_Low', 'Obs_Low_High', 'Obs_High_High'
        #         , 'ChiPval', 'Chi2'
    ])
    return test_res, random_test_res

def stat_test_plot(test_res, random_test_res, x_max = 100, y_max = 1):
    plt.figure(figsize=(10,5))
    for i, Obs in enumerate(['Obs_Low_Low', 'Obs_High_Low', 'Obs_Low_High', 'Obs_High_High']):
        plt.subplot(221 + i)
        plt.hist([
            np.array(test_res[Obs], dtype = int),
            np.array(random_test_res[Obs], dtype = int),],
            label=['SL','Random'], histtype=u'step',  density = True)
        plt.legend(loc='upper right')
        plt.xlabel(Obs)
        plt.ylabel('Frequency')
        plt.xlim(-1, x_max)
        plt.ylim(0, y_max)

    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.5,wspace=0.35)

    plt.show()
