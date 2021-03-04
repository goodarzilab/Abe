import pandas as pd
from glob import glob

def read_gmt(PATH):
    '''
    Read given gmt file into a Dictionary 
    '''
    with open(PATH) as gmt:
        lines = gmt.readlines()
        out = {}
        for line in lines:
            data = line.split('\t')
            name = data[0]
            url  = data[1]
            genes= data[2:]
            genes[-1] = genes[-1].split('\n')[0]

            out[name] = {}
            out[name]['url'] = url 
            out[name]['genes'] = genes
            
    return out


def read_page_index(PATH):
    '''
    Read given *index.txt files into a Dictionary 
    '''
    with open(PATH) as raw:
        lines = raw.readlines()
        genes_dict = {}
        for line in lines:
            data = line.split('\t')
            gene = data[0]
            genesets  = data[1:]
            genesets[-1] = genesets[-1].split('\n')[0]

            genes_dict[gene] = genesets
            
    # list all unique genesets 
    all_gs = []
    for gs in genes_dict.values():
        all_gs = all_gs + gs
    all_gs = set(all_gs)
    
    # dictionary which keys are genesets and values are associated genes 
    gs_dict = {}
    for gs in all_gs:
        gs_dict[gs] = {g for g in genes_dict.keys() if gs in genes_dict[g]}
        
    out = gs_dict
    
    return out


def read_page_names(PATH):
    '''
    Read given *name.txt files into a Dictionary 
    '''
    with open(PATH) as raw:
        lines = raw.readlines()
        out = {}
    for line in lines:
            data = line.split('\t')
            name0= data[0]
            name1= data[1]
            pw_type= data[2].split('\n')[0]
            out[name0] = [name1, pw_type]
    return out


def read_page_annotations(gs_name,ANNDIR='/flash/bin/iPAGEv1.0/PAGE_DATA/ANNOTATIONS/'):
    '''
    Read gene set annotations into python from PAGE_DATA format
    '''
    annotations = {}
    gmt = glob(f'{ANNDIR}{gs_name}/*.gmt')
    index = glob(f'{ANNDIR}{gs_name}/*_index.txt')
    names = glob(f'{ANNDIR}{gs_name}/*_names.txt')
    if index:
        annotations['index'] = read_page_index(index[0])
#     if names:
#         annotations['names'] = read_page_names(names[0])
#     if gmt:
#         annotations['gmt'] = read_gmt(gmt[0])

    return annotations


def make_page_dict(PATH):
    '''
    PATH = a complete path to a pvmatrix.txt file, part of results from iPAGE run 
    
    Processes: 
    1) Read p-value matrix data into a data frame
    2) Include annotations for the gene set from the PAGE directory
    
    Output: Python dictionary contain pvmatrix and related annotations to the gene set
    '''
    ### 1 ### 
    # read pvmatrix.txt file 
    df = pd.read_csv(PATH, sep='\t',index_col=0)
    # remove duplicated named (row) names 
    if all([geneset.split(' ')[0] == geneset.split(' ')[1] for geneset in df.index.tolist()]):
        df.index = [geneset.split(' ')[0] for geneset in df.index.tolist() ]
        
    ### 2 ### 
    gs_name = PATH.split('/')[-2]
    ann = read_page_annotations(gs_name)

    out = {}
    out['gs_name'] = gs_name
    out['annotations'] = ann
    out['data'] = df
    
    return out
