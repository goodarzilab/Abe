import os
import ntpath
import subprocess
import glob
import argparse
import warnings
import re
import pandas as pd
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")

cnt_ex = glob.glob('/rumi/shams/abe/People/Alex_Ge/bam/*_exons.txt')[0]
cnt_in = glob.glob('/rumi/shams/abe/People/Alex_Ge/bam/*_introns.txt')[0]
tmp_ex = pd.read_table(cnt_ex, skiprows=1)['Geneid'].tolist()
tmp_in = pd.read_table(cnt_in, skiprows=1)['Geneid'].tolist()
df_ex = pd.DataFrame({}, index=tmp_ex)
df_in = pd.DataFrame({}, index=tmp_in)

files_ex = glob.glob(ntpath.dirname(cnt_ex)+'/*_exons.txt')
files_in = glob.glob(ntpath.dirname(cnt_in)+'/*_introns.txt')
for f in files_ex:
    counts = pd.read_table(f, skiprows=1).iloc[:, -1].tolist()
    f_n = f.split('/')[-1].replace('_merge.bam_exons.txt', '')
    df_ex[f_n] = counts
for f in files_in:
    counts = pd.read_table(f, skiprows=1).iloc[:, -1].tolist()
    f_n = f.split('/')[-1].replace('_merge.bam_introns.txt', '')
    df_in[f_n] = counts
df_ex = df_ex[df_in.columns]


print('\n************************* Generate REMBRANDTS Inputs *************************\n')

# Prepare REMBRANDTS inputs
stbldir = '/rumi/shams/abe/People/Alex_Ge/stbl/'
if not (os.path.exists(stbldir)):
    os.mkdir(stbldir)

ol = list(set(df_ex.index)&set(df_in.index))
df_ex = df_ex.loc[ol, ]
df_in = df_in.loc[ol, ]
#
countdir = stbldir+'count_tables/'
if not (os.path.exists(countdir)):
    os.mkdir(countdir)
for c in df_ex.columns:
    df_ex.loc[:, c].to_csv(countdir+'%s.exonic.txt'%c, sep='\t', header=None)
    df_in.loc[:, c].to_csv(countdir+'%s.intronic.txt'%c, sep='\t', header=None)

labels = np.repeat(df_ex.columns.tolist(), 2)
files = [labels[i]+'.exonic.txt' if i%2==0 else labels[i]+'.intronic.txt' for i in range(len(labels))]
readtypes = ['exonic' if i%2==0 else 'intronic' for i in range(len(labels))]
meta = pd.DataFrame({'Label':labels, 'File':files, 'ReadType':readtypes, 'Batch':1})[['Label', 'File', 'ReadType', 'Batch']]
meta.to_csv(stbldir+'REMBRANDTS_input.txt', index=False, sep='\t')

print('\n******************************* Run REMBRANDTS *******************************\n')

# Run REMBRANDTS
os.chdir('/flash/hani/bin/REMBRANDTS')
subprocess.call('bash REMBRANDTS.sh %s %s %s 0.99 linear'%('Alex_Ge', stbldir+'REMBRANDTS_input.txt', countdir), shell=True)
if not os.path.exists(stbldir+'out/'):
  os.mkdir(stbldir+'out/')
subprocess.call('mv ./out/%s/ %sout/' % ('Alex_Ge', stbldir), shell=True)
# subprocess.call('rm -r ./out/%s' % ('Alex_Ge'), shell=True)
subprocess.call('mv ./tmp/%s/ %stmp/' % ('Alex_Ge', stbldir), shell=True)
# subprocess.call('rm -r ' + REMBRANDTS_dir +'/tmp/' % ('Alex_Ge'), shell=True)

out = pd.read_table(stbldir+'out/%s/stability.filtered.mx.txt'%'Alex_Ge', index_col=0)
out.columns = [i.replace('.x','').replace('.','-') for i in out.columns]
if 'X' in out.columns[0] and 'X' not in df_ex.columns[0]:
  out.columns = [i[1:] for i in out.columns]
out.to_csv(stbldir+'REMBRANDTS_stbl.txt', sep='\t')

print('\n********************************** Finished **********************************\n')
