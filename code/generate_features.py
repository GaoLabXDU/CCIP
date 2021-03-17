#!/usr/bin/env python
# coding: utf-8

import os,sys
import numpy as np
import pandas as pd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os,sys
import networkx as nx
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-o', '--output_path',default='./output/',help='Path for output')
parser.add_option('-c', '--ctcf_file',default='./CTCF_peak.bed',help='CTCF ChIP-seq data')
parser.add_option('-m', '--ctcf_motif_file',default='./fimo.csv',help='CTCF motif occurence data')
parser.add_option('-p', '--chia_pet_file',default='ctcf.interactions.intra.bedpe',help='CTCF ChIA-PET data')
parser.add_option('-r', '--rad21_file',default='./rad21.narrowPeak',help='RAD21 ChIP-seq data')
parser.add_option('-a', '--age_file',default='./CTCF_age.bed',help='CTCF age data')
(opts, args) = parser.parse_args()

output_path=opts.output_path
if not output_path.endswith('/'):
    output_path=output_path+'/'
if not os.path.exists(output_path):
    os.makedirs(output_path)

chia_pet_file=opts.chia_pet_file

ctcf_file=opts.ctcf_file

ctcf_motif_file=opts.ctcf_motif_file


rad21_file=opts.rad21_file

age_file=opts.age_file


# # motif data
motif=pd.read_csv(ctcf_motif_file)
motif=motif.rename(columns={'chromosome':'chr','sth1':'score','sth2':'pValue','sth3':'qValue'})
motif["name"]=motif['name']=motif['chr']+':'+motif['start'].astype(str)+'-'+motif['end'].astype(str)
motif.columns=['motif_'+name for name in list(motif.columns)]
motif.motif_strand=motif.motif_strand.transform(lambda x:1 if x=='+' else 0)
motif[['motif_chr','motif_start','motif_end','motif_name']].to_csv(output_path+'ctcf_motif.bed','\t',header=False,index=False)


# # CTCF data
ctcf=pd.read_csv(ctcf_file,sep='\t',header=None,names=['chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'])
ctcf['name']=ctcf['chr']+':'+ctcf['start'].astype(str)+'-'+ctcf['end'].astype(str)
ctcf.peak=ctcf.start+ctcf.peak
ctcf.columns=['ctcf_'+name for name in list(ctcf.columns)]
ctcf[['ctcf_chr','ctcf_start','ctcf_end','ctcf_name']].to_csv(output_path+'ctcf.bed',sep='\t',header=False,index=False)

os.system('bedtools intersect -f 1 -a %sctcf_motif.bed -b %sctcf.bed -wa -wb>%smotif_to_ctcf_intersection.bed'%(output_path,output_path,output_path))

motif_to_ctcf=pd.read_csv(output_path+'motif_to_ctcf_intersection.bed','\t',header=None,usecols=[3,7],names=['motif_name','ctcf_name'])

anchor=motif.merge(motif_to_ctcf).merge(ctcf)

anchor=anchor.drop_duplicates(subset=['motif_name'],keep='first')


# # rad21

rad=pd.read_csv(rad21_file,sep='\t',header=None,names=['chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'])
rad['name']=rad['chr']+':'+rad['start'].astype(str)+'-'+rad['end'].astype(str)
rad.peak=rad.start+rad.peak
rad.columns=['rad_'+name for name in list(rad.columns)]
rad[['rad_chr','rad_start','rad_end','rad_name']].to_csv(output_path+'rad21.bed',sep='\t',header=False,index=False)

os.system('bedtools intersect -f 0.1 -a %sctcf_motif.bed -b %srad21.bed -wa -wb>%smotif_to_rad21_intersection.bed'%(output_path,output_path,output_path))

motif_to_rad=pd.read_csv(output_path+'motif_to_rad21_intersection.bed','\t',header=None,usecols=[3,7],names=['motif_name','rad_name'])

anchor=anchor.merge(motif_to_rad,how='left').merge(rad,how='left')
anchor=anchor.drop_duplicates(subset=['motif_name'],keep='first')

age=pd.read_csv(age_file,sep='\t')
age=age.rename(columns={'chromosome':'chr'})
age['name']=age['chr']+':'+age['start'].astype(str)+'-'+age['end'].astype(str)
age.columns=['age_'+name for name in list(age.columns)]
age=age.rename(columns={'age_name':'motif_name'})

anchor=anchor.merge(age)

def age_map(age_name):
    age_dict={'hg19':1,'Human-Chimp':2,'Homininae':3,'Hominidae':4,
	'Catarrhini':5,'Simiiformes':6,'Haplorrhini':7,'Primate':8,
	'Strepsirrhini':9,'Human-Mouse':10,'Node_11':11,'Node_12':12,
	'Node_13':13,'Node_14':14,'Node_15':15,'Node_16':16,'Node_17':17,'Root':18,'NA':0,np.nan:0}
    return age_dict[age_name]


anchor['age']=anchor.age_age.transform(age_map)

anchor=anchor[['motif_name','motif_chr','motif_start','motif_end','motif_strand','motif_score','motif_motif','ctcf_signalValue','rad_signalValue','age']]

anchor=anchor.fillna(0)

anchor=anchor.sort_values(by=['motif_chr','motif_start'],ascending=True)

anchor['motif_index']=anchor.index.astype('object')

anchor['distance_last']=0
anchor['distance_next']=0
for chrom in anchor.motif_chr.unique():
    anchor_part=anchor[anchor.motif_chr==chrom]
    anchor.loc[anchor_part.index[:-1],'distance_last']=anchor_part.motif_start.values[1:]-anchor_part.motif_start.values[:-1]
    anchor.loc[anchor_part.index[1:],'distance_next']=anchor_part.motif_start.values[1:]-anchor_part.motif_start.values[:-1]

def one_hot(row):
    seq=row.motif_motif
    seq_dict = {'A':[1, 0, 0, 0], 'G':[0, 1, 0, 0],
                'C':[0, 0, 1, 0], 'T':[0, 0, 0, 1],
                'a':[1, 0, 0, 0], 'g':[0, 1, 0, 0],
                'c':[0, 0, 1, 0], 't':[0, 0, 0, 1]}
    temp = []
    for c in seq:
        temp.extend(seq_dict.get(c, [0, 0, 0, 0]))
    return temp

column_name=[c+str(num) for num in range(19) for c in ['A','G','C','T']]

anchor_motif=pd.DataFrame(np.array(list(anchor.apply(one_hot,axis=1).values)),columns=column_name,index=anchor.index)

anchor=pd.concat([anchor,anchor_motif],axis=1)

anchor.to_csv(output_path+'anchor.csv',index=False)

pair=pd.read_csv(output_path+'pair_all_balance.csv')

from_anchor_columns=['from_'+name for name in  anchor.columns]
to_anchor_columns=['to_'+name for name in  anchor.columns]
from_anchor=anchor.copy()
to_anchor=anchor.copy()
from_anchor.columns=from_anchor_columns
to_anchor.columns=to_anchor_columns

pair=pair.merge(from_anchor)
pair=pair.merge(to_anchor)


def get_positive_motif_apply(row):
    between_motif=anchor.loc[row.from_motif_index+1:row.to_motif_index-1]
    return np.sum(between_motif.motif_strand==1)
def get_negative_motif_apply(row):
    between_motif=anchor.loc[row.from_motif_index+1:row.to_motif_index-1]
    return np.sum(between_motif.motif_strand==0)
def get_motif_score_apply(row):
    between_motif=anchor.loc[row.from_motif_index+1:row.to_motif_index-1]
    return np.sum(between_motif.motif_score)
def get_positive_motif_score_apply(row):
    between_motif=anchor.loc[row.from_motif_index+1:row.to_motif_index-1]
    return np.sum(between_motif[between_motif.motif_strand==1].motif_score)
def get_negative_motif_score_apply(row):
    between_motif=anchor.loc[row.from_motif_index+1:row.to_motif_index-1]
    return np.sum(between_motif[between_motif.motif_strand==0].motif_score)

def get_ctcf_signal_apply(row):
    between_motif=anchor.loc[row.from_motif_index+1:row.to_motif_index-1]
    return np.sum(between_motif.ctcf_signalValue)
def get_positive_ctcf_signal_apply(row):
    between_motif=anchor.loc[row.from_motif_index+1:row.to_motif_index-1]
    return np.sum(between_motif[between_motif.motif_strand==1].ctcf_signalValue)
def get_negative_ctcf_signal_apply(row):
    between_motif=anchor.loc[row.from_motif_index+1:row.to_motif_index-1]
    return np.sum(between_motif[between_motif.motif_strand==0].ctcf_signalValue)


pair['between_motif']=pair.to_motif_index-pair.from_motif_index
pair['between_positive_motif']=pair.apply(get_positive_motif_apply,axis=1)
pair['between_negative_motif']=pair.apply(get_negative_motif_apply,axis=1)

pair['between_ctcf_score']=pair.apply(get_motif_score_apply,axis=1)
pair['between_positive_ctcf_score']=pair.apply(get_positive_motif_score_apply,axis=1)
pair['between_negative_ctcf_score']=pair.apply(get_negative_motif_score_apply,axis=1)

pair['between_ctcf_signalValue']=pair.apply(get_ctcf_signal_apply,axis=1)
pair['between_positive_ctcf_signalValue']=pair.apply(get_positive_ctcf_signal_apply,axis=1)
pair['between_negative_ctcf_signalValue']=pair.apply(get_negative_ctcf_signal_apply,axis=1)


pair.to_csv(output_path+'sample.csv',index=False)