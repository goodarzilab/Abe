import numpy as np
import pandas as pd 
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import scipy as sp

#adapted from https://github.com/mhorlbeck/ScreenProcessing
figureScale = 2

def cleanAxes(axis, top=False, right=False, bottom=True, left=True):
    axis.spines['top'].set_visible(top)
    axis.spines['right'].set_visible(right)
    axis.spines['left'].set_visible(left)
    axis.spines['bottom'].set_visible(bottom)

    #turn off all ticks
    axis.yaxis.set_ticks_position('none')
    axis.xaxis.set_ticks_position('none')

    #now re-enable visibles
    if top:
        axis.xaxis.tick_top()
    if bottom:
        axis.xaxis.tick_bottom()
    if left:
        axis.yaxis.tick_left()
    if right:
        axis.yaxis.tick_right()

def plotGrid(axis, vert_origin = True, horiz_origin=True, unity=True):
    ylim = axis.get_ylim()
    xlim = axis.get_xlim()
    if vert_origin:
        axis.plot((0,0), ylim, color='#BFBFBF', lw=.5, alpha=.5) 
    if horiz_origin:
        axis.plot(xlim,(0,0), color='#BFBFBF', lw=.5, alpha=.5) 
    if unity:
        xmin = min(xlim[0], ylim[0])
        xmax = max(xlim[1], ylim[1])
        axis.plot((xmin,xmax),(xmin,xmax), color='#BFBFBF', lw=.5, alpha=.5) 
        
    axis.set_ylim(ylim)
    axis.set_xlim(xlim)

def displayFigure(fig, savetitle=None,imageExtension='pdf'):
    if savetitle:
        plt.show(fig)
        
        fullTitle =  os.path.join('{0}.{1}'.format(savetitle, imageExtension))
        print(fullTitle)
        fig.savefig(fullTitle, dpi=1000)
        plt.close(fig) 
        return fullTitle        
    else:
        plt.show(fig)
        
##gene-level plotting functions
def volcanoPlot(table, effectSizeLabel='log2FoldChange', pvalueLabel='pvalue', hitThreshold=7, labelHits = True):

    discScore = lambda z,p: p * np.abs(z)
    
    table.loc[:,'thresh'] = discScore(table[effectSizeLabel],-1*np.log10(table[pvalueLabel])) >= hitThreshold

    yGenes = -1*np.log10(table[pvalueLabel])
    xGenes = table[effectSizeLabel]

    fig, axis = plt.subplots(1,1, figsize=(4*figureScale,3.5*figureScale))
    cleanAxes(axis)

    axis.scatter(table.loc[table['thresh'],effectSizeLabel], -1*np.log10(table.loc[table['thresh'],pvalueLabel].values),
                 s=4,
                 c='#7570b3',
                 label = 'Gene hit',
                 rasterized=True)
                 
    if labelHits:
        for gene, row in table.loc[table['thresh']].iterrows():
            axis.text(row[effectSizeLabel], -1*np.log10(row[pvalueLabel]), gene, fontsize=6,
            horizontalalignment = 'left' if row[effectSizeLabel] > 0 else 'right', verticalalignment='center')
            

    plotGrid(axis, vert_origin=True, horiz_origin=False, unity=False)

    ymax = np.ceil(max(yGenes)) * 1.02
    xmin = min(xGenes) * 1.05
    xmax = max(xGenes) * 1.05
    
    axis.plot(np.linspace(xmin,xmax,1000),np.abs(hitThreshold/np.linspace(xmin,xmax,1000)),'k--', lw=.5)

    axis.set_xlim((xmin,xmax))
    axis.set_ylim((0,ymax))

#     axis.set_xlabel('{3} {0} {1} ({2})'.format(phenotype, replicate, effectSizeLabel, 'gene' if not transcripts else 'transcript'),fontsize=8)
    axis.set_xlabel('{0}'.format(effectSizeLabel, 'gene' ),fontsize=8)
    axis.set_ylabel('-log10 {0}'.format(pvalueLabel),fontsize=8)
    
    plt.legend(loc='best', fontsize=6, handletextpad=0.005)

    plt.tight_layout()
    return displayFigure(fig)
