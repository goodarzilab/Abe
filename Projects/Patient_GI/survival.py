import pandas as pd
import numpy as np
from time import time
# import lifelines as ll


# Read combined study clinical data (downloaded from cBioPortal):
def SL_survival(data, test):
    t0 = time()
    clinical_data = pd.read_csv('cBioPortal/combined_study_clinical_data.csv')
    dfs = []
    for t in test:
        survival_data = []
        g1, g2 = t.split('_')
        g1_index = [i for i, g in enumerate(data['raw']['Genes']) if g == g1][0]
        g2_index = [i for i, g in enumerate(data['raw']['Genes']) if g == g2][0]
        Q1g1 = np.quantile(data['raw']['RNA-Seq'][g1_index], 0.25)  # Gene 1 threshold
        Q1g2 = np.quantile(data['raw']['RNA-Seq'][g2_index], 0.25)  # Gene 2 threshold
        for p, sam in enumerate(data['raw']['sample_ids']):
            dat1 = 1 * (data['raw']['RNA-Seq'][data['raw']['Genes'] == g1][p] <= Q1g1)
            dat2 = 1 * (data['raw']['RNA-Seq'][data['raw']['Genes'] == g2][p] <= Q1g2)
            surS = clinical_data['Overall Survival Status'][clinical_data['Sample ID'] == sam]
            surM = clinical_data['Overall Survival (Months)'][clinical_data['Sample ID'] == sam]
            if len(surS) > 0:
                survival_data.append([sam, dat1, dat2, surS.values[0], surM.values[0]])
        dfs.append(
            pd.DataFrame(survival_data, columns=['sample_ids', g1 + ' is low', g2 + ' is low',
                                                 'Status', 'Months']).dropna()
        )
    print("done in %fs" % (time() - t0))
    return dfs
