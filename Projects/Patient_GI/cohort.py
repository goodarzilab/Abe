import pandas as pd
import numpy as np
from time import time


# make final data set
def make_study_dataset(f_path, get_RNASeq=False):
    t0 = time()
    with open(f_path) as fp:
        lines = [l.split('\t') for l in fp.readlines()]
        # switch NAs -> 0s
        for i, l in enumerate(lines):
            if 'NA' in l:
                n = len(l[2:])
                l = l[0:2]
                for z in [0] * n:
                    l.append(z)
            lines[i] = l
        # genes = [d[0] for d in lines[1:]]
        RNA_seq = [np.array(d[2:], dtype=float) for d in lines[1:]]

        Q1 = [np.quantile(p, 0.1) for p in RNA_seq]  # Low expression threshold
        Q3 = [np.quantile(p, 0.9) for p in RNA_seq]  # High expression threshold

    data = pd.DataFrame(data=[i + j for i, j in zip([-1 * (p <= q1) for p, q1 in zip(RNA_seq, Q1)],
                                                    [1 * (p >= q3) for p, q3 in zip(RNA_seq, Q3)])],
                        index=[d[0] for d in lines[1:]], columns=lines[0][2:])

    print("done in %fs" % (time() - t0))
    if not get_RNASeq:
        return data
    else:
        RNA_seq = pd.DataFrame(data=RNA_seq, index=[d[0] for d in lines[1:]],
                               columns=lines[0][2:])
        return RNA_seq, data
