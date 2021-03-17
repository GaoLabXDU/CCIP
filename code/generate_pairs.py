#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os,sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-o', '--output_path',default='./output/',help='Path for output')
parser.add_option('-c', '--ctcf_file',default='./CTCF_peak.bed',help='CTCF ChIP-seq data')
parser.add_option('-m', '--ctcf_motif_file',default='./fimo.csv',help='CTCF motif occurence data')
parser.add_option('-p', '--chia_pet_file',default='ctcf.interactions.intra.bedpe',help='CTCF ChIA-PET data')
(opts, args) = parser.parse_args()

output_path=opts.output_path
if not output_path.endswith('/'):
    output_path=output_path+'/'
if not os.path.exists(output_path):
    os.makedirs(output_path)

chia_pet_file=opts.chia_pet_file

ctcf_file=opts.ctcf_file

ctcf_motif_file=opts.ctcf_motif_file


#ChIA-PET data

gm12878=pd.read_csv(chia_pet_file,header=None,sep='\t',names=['chr','from_start_cord','from_end_cord'
                                                                  ,'to_chr','to_start_cord','to_end_cord','from_name','to_name','from_score','to_score'
                                                                  ,'frequency'],usecols=['chr','from_start_cord','from_end_cord'
                                                                  ,'to_start_cord','to_end_cord'
                                                                  ,'frequency'])


gm12878['from_length']=gm12878['from_end_cord']-gm12878['from_start_cord']
gm12878['to_length']=gm12878['to_end_cord']-gm12878['to_start_cord']
gm12878['span']=(gm12878['to_start_cord']+gm12878['to_end_cord'])//2-(gm12878['from_start_cord']+gm12878['from_end_cord'])//2
gm12878['window_size']=gm12878['to_start_cord']-gm12878['from_end_cord']
gm12878['ratio']=gm12878['window_size']/(gm12878['from_length']+gm12878['to_length'])

gm12878['from_name']=gm12878['chr']+':'+gm12878['from_start_cord'].astype(str)+'-'+gm12878['from_end_cord'].astype(str)
gm12878['to_name']=gm12878['chr']+':'+gm12878['to_start_cord'].astype(str)+'-'+gm12878['to_end_cord'].astype(str)

# ## CTCF data
ctcf=pd.read_csv(ctcf_file,sep='\t',header=None,names=['chr','start','end','name','score','strand','signalValue','pValue','qValue','peak'])
ctcf['name']=ctcf['chr']+':'+ctcf['start'].astype(str)+'-'+ctcf['end'].astype(str)
ctcf[['chr','start','end','name']].to_csv(output_path+'ctcf.bed',sep='\t',header=False,index=False)

# ## CTCF motif data

motif_ctcf=pd.read_csv(ctcf_motif_file)
motif_ctcf=motif_ctcf.drop(columns=['name'])
motif_ctcf=motif_ctcf.rename(columns={'chromosome':'chr'})
motif_ctcf["name"]=motif_ctcf['name']=motif_ctcf['chr']+':'+motif_ctcf['start'].astype(str)+'-'+motif_ctcf['end'].astype(str)

motif_ctcf[['chr','start','end','name']].to_csv(output_path+'ctcf_motif.bed',sep='\t',header=None,index=False)


# # Anchor list of ChIA-PET data

anchors_left=gm12878[['chr','from_start_cord','from_end_cord']].copy()
anchors_left.columns=['chr','start','end']
anchors_right=gm12878[['chr','to_start_cord','to_end_cord']].copy()
anchors_right.columns=['chr','start','end']
anchors=pd.concat([anchors_left,anchors_right],axis=0)

anchors.sort_values(by=['chr','start','end'],inplace=True)
anchors['name']=anchors['chr']+':'+anchors['start'].astype(str)+'-'+anchors['end'].astype(str)

anchors=anchors.drop_duplicates()

anchors.to_csv(output_path+'ctcf_pet_anchors.bed',sep='\t',header=False,index=False)


# # Interactions between different data

os.system('bedtools intersect -F 0.1 -a %sctcf_pet_anchors.bed -b %sctcf.bed -wa -wb >%sanchor_to_ctcf_intersection.bed'%(output_path,output_path,output_path))
os.system('bedtools intersect -F 0.1 -a %sctcf.bed -b %sctcf_motif.bed -wa -wb>%sctcf_to_motif_intersection.bed'%(output_path,output_path,output_path))
#os.system('bedtools intersect -F 0.1 -a %sctcf.bed -b %sctcf_pet_anchors.bed -wa -wb >%sctcf_to_anchor_intersection.bed'%(output_path,output_path,output_path))

anchor_to_ctcf=pd.read_csv(output_path+'anchor_to_ctcf_intersection.bed',header=None,sep='\t',names=['anchor_name','ctcf_name'],usecols=[3,7])

ctcf_to_motif=pd.read_csv(output_path+'ctcf_to_motif_intersection.bed','\t',header=None,usecols=[3,4,5,6,7],names=['ctcf_name','motif_chr','motif_start','motif_end','motif_name'])


# # generate positive samples

gm12878_extend=(gm12878
                .merge(anchor_to_ctcf,left_on='from_name',right_on='anchor_name',how='inner')
                .rename({'ctcf_name':'from_ctcf_name'},axis=1)
                .drop('anchor_name',axis=1)
                .merge(anchor_to_ctcf,left_on='to_name',right_on='anchor_name',how='inner')
                .rename({'ctcf_name':'to_ctcf_name'},axis=1)
                .drop('anchor_name',axis=1)
                .merge(ctcf_to_motif,left_on='from_ctcf_name',right_on='ctcf_name',how='inner')
                .rename({'motif_name':'from_motif_name','motif_chr':'from_motif_chr','motif_start':'from_motif_start','motif_end':'from_motif_end'},axis=1)
                .drop('ctcf_name',axis=1)
                .merge(ctcf_to_motif,left_on='to_ctcf_name',right_on='ctcf_name',how='inner')
                .rename({'motif_name':'to_motif_name','motif_chr':'to_motif_chr','motif_start':'to_motif_start','motif_end':'to_motif_end'},axis=1)
                .drop('ctcf_name',axis=1)
               )
gm12878_extend['ctcf_pair_name']=gm12878_extend['from_ctcf_name'].astype(str)+'+'+gm12878_extend['to_ctcf_name'].astype(str)

gm12878_one2one=gm12878_extend.drop_duplicates(subset=['from_name','to_name'],keep=False).copy()

gm12878_m2m=gm12878_extend[gm12878_extend.duplicated(subset=['from_name','to_name'],keep=False)].copy()

positive_pair_set=set(gm12878_one2one['ctcf_pair_name'])|set(gm12878_m2m['ctcf_pair_name'])


# ## filter out samples by distance

pair_positive=gm12878_one2one[['chr','from_name','to_name','from_ctcf_name','to_ctcf_name','from_motif_name','from_motif_chr','from_motif_start','from_motif_end','to_motif_name','to_motif_chr','to_motif_start','to_motif_end','frequency']].copy()
pair_positive=pair_positive.merge(motif_ctcf[['name','start','end']],left_on='from_motif_name',right_on='name').drop('name',axis=1)
pair_positive=pair_positive.merge(motif_ctcf[['name','start','end']],left_on='to_motif_name',right_on='name').drop('name',axis=1)
pair_positive['from_cord']=pair_positive.eval('(end_x+start_x)/2')
pair_positive['to_cord']=pair_positive.eval('(end_y+start_y)/2')
pair_positive=pair_positive.drop(['start_x','start_y','end_x','end_y'],axis=1).copy()

pair_positive['from_cord']=pair_positive['from_cord'].astype(np.int32)
pair_positive['to_cord']=pair_positive['to_cord'].astype(np.int32)
pair_positive['distance']=pair_positive.eval('to_cord-from_cord')
pair_positive=pair_positive[pair_positive.eval('10000<=distance<=1000000')]


# # generate negative samples

ctcf_extend=ctcf.merge(ctcf_to_motif,left_on='name',right_on='ctcf_name',how='inner').drop('ctcf_name',axis=1)

ctcf_one2one=ctcf_extend.drop_duplicates(subset=['name'],keep=False)

ctcf_one2one=ctcf_one2one.drop_duplicates(subset=['motif_name'],keep=False)

pair_candidate=ctcf_one2one.merge(ctcf_one2one,left_on='chr',right_on='chr')

pair_candidate['ctcf_pair_name']=pair_candidate['name_x'].astype(str)+'+'+pair_candidate['name_y'].astype(str)

index_filt_out=pair_candidate['ctcf_pair_name'].apply(lambda x:x in positive_pair_set)
pair_candidate_valid=pair_candidate[np.logical_not(index_filt_out)].copy()
pair_negative=pair_candidate_valid[['chr','motif_name_x','motif_name_y','motif_chr_x','motif_start_x','motif_end_x','motif_chr_y','motif_start_y','motif_end_y']].copy()

pair_negative['from_cord']=pair_candidate_valid.eval('(motif_end_x+motif_start_x)/2')
pair_negative['to_cord']=pair_candidate_valid.eval('(motif_end_y+motif_start_y)/2')

pair_negative['from_cord']=pair_negative['from_cord'].astype(np.int32)
pair_negative['to_cord']=pair_negative['to_cord'].astype(np.int32)
pair_negative['distance']=pair_negative.eval('to_cord-from_cord')

pair_negative=pair_negative[pair_negative.eval('10000 <= distance <=1000000')].copy()

pair_negative.rename({'motif_name_x':'from_motif_name','motif_chr_x':'from_motif_chr','motif_start_x':'from_motif_start','motif_end_x':'from_motif_end','motif_name_y':'to_motif_name','motif_chr_y':'to_motif_chr','motif_start_y':'to_motif_start','motif_end_y':'to_motif_end'},axis=1,inplace=True)

pair_negative=pair_negative.merge(ctcf_extend[['name','motif_name']],left_on='from_motif_name',right_on='motif_name').drop('motif_name',axis=1).rename({'name':'from_ctcf_name'},axis=1)
pair_negative=pair_negative.merge(ctcf_extend[['name','motif_name']],left_on='to_motif_name',right_on='motif_name').drop('motif_name',axis=1).rename({'name':'to_ctcf_name'},axis=1)

pair_negative['from_name']='null'
pair_negative['to_name']='null'
pair_negative['frequency']=0
pair_negative=pair_negative[pair_positive.columns]


pair_positive.reset_index(drop=True,inplace=True)
pair_negative.reset_index(drop=True,inplace=True)

pair_positive['bin'],bins=pd.qcut(pair_positive['distance'],10,precision=1,retbins=True,labels=range(10))
pair_negative['bin']=pd.cut(pair_negative['distance'],bins,precision=1,retbins=False,include_lowest=True,labels=range(10))


# # extract motif features
pair_positive=(pair_positive
               .merge(motif_ctcf[['name','strand']],left_on='from_motif_name',right_on='name',how='left')
               .drop('name',axis=1)
               .rename({'strand':'from_ctcf_motif_strand'},axis=1)
               .merge(motif_ctcf[['name','strand']],left_on='to_motif_name',right_on='name',how='left')
               .drop('name',axis=1)
               .rename({'strand':'to_ctcf_motif_strand'},axis=1)
              )
pair_negative=(pair_negative
               .merge(motif_ctcf[['name','strand']],left_on='from_motif_name',right_on='name',how='left')
               .drop('name',axis=1)
               .rename({'strand':'from_ctcf_motif_strand'},axis=1)
               .merge(motif_ctcf[['name','strand']],left_on='to_motif_name',right_on='name',how='left')
               .drop('name',axis=1)
               .rename({'strand':'to_ctcf_motif_strand'},axis=1)
              )

pair_positive_convergence=pair_positive[(pair_positive.from_ctcf_motif_strand=='+')&(pair_positive.to_ctcf_motif_strand=='-')]
pair_negative_convergence=pair_negative[(pair_negative.from_ctcf_motif_strand=='+')&(pair_negative.to_ctcf_motif_strand=='-')]
pair_positive_tandem=pair_positive[pair_positive.from_ctcf_motif_strand==pair_positive.to_ctcf_motif_strand]
pair_negative_tandem=pair_negative[pair_negative.from_ctcf_motif_strand==pair_negative.to_ctcf_motif_strand]

positive=pair_positive_tandem.copy()
negative=pair_negative_tandem.copy()

frequency_threshold=2
num=10
positive=positive[positive.frequency>=frequency_threshold]
positive['bin'],bins=pd.qcut(positive['distance'],num,precision=1,retbins=True,labels=range(num))
negative['bin']=pd.cut(negative['distance'],bins,precision=1,retbins=False,include_lowest=True,labels=range(num))

positive_bin_size=positive[['frequency','bin']].groupby('bin').count()
negative_bin_size=negative[['frequency','bin']].groupby('bin').count()

neg_sample_list=[]
for name,group in negative.groupby(negative.bin,as_index=False):
    p_num=positive_bin_size.iloc[name,0]
    n_num=negative_bin_size.iloc[name,0]
    neg_sample_list.append(group.sample(p_num,replace=False if p_num<=n_num else True).reset_index(drop=True))
neg_sample=pd.concat(neg_sample_list)
sample_tandem=pd.concat([positive,neg_sample])

#sample=sample.drop(['bin'],axis=1)
sample_tandem=sample_tandem.reset_index(drop=True)
sample_tandem=sample_tandem.sample(frac=1,random_state=0).reset_index(drop=True)
#sample_tandem.to_csv(output_path+'pair_tandem.csv',index=False)

positive=pair_positive_convergence.copy()
negative=pair_negative_convergence.copy()

frequency_threshold=2
num=10
positive=positive[positive.frequency>=frequency_threshold]
positive['bin'],bins=pd.qcut(positive['distance'],num,precision=1,retbins=True,labels=range(num))
negative['bin']=pd.cut(negative['distance'],bins,precision=1,retbins=False,include_lowest=True,labels=range(num))

positive_bin_size=positive[['frequency','bin']].groupby('bin').count()
negative_bin_size=negative[['frequency','bin']].groupby('bin').count()

neg_sample_list=[]
for name,group in negative.groupby(negative.bin,as_index=False):
    p_num=positive_bin_size.iloc[name,0]
    n_num=negative_bin_size.iloc[name,0]
    neg_sample_list.append(group.sample(p_num,replace=False if p_num<=n_num else True).reset_index(drop=True))
neg_sample=pd.concat(neg_sample_list)
sample_convergence=pd.concat([positive,neg_sample])

#sample=sample.drop(['bin'],axis=1)
sample_convergence=sample_convergence.reset_index(drop=True)
sample_convergence=sample_convergence.sample(frac=1,random_state=0).reset_index(drop=True)
#sample_convergence.to_csv(output_path+'pair_convergence.csv',index=False)

sample=pd.concat([sample_convergence,sample_tandem]).reset_index(drop=True)
sample.to_csv(output_path+'pair_all_balance.csv',index=False)
