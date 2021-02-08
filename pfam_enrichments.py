# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 13:26:17 2020

@author: asafl-lab
"""

import pandas as pd
import glob
import scipy.stats as stats
import numpy as np
from statsmodels.stats.multitest import multipletests



names = pd.read_table('all_names.tsv')

ecis = pd.read_table('newest_db_tabbed.tsv')
ecis.gene_oid = ecis.gene_oid.astype('int64')

gene_in_ecis = ecis.gene_oid.to_list()

column_names = ['pfam_id', 'pfam_name','num_pfams_in_eCIS', 'all_pfams_counts','num_pfam_NOT_in_ecis', 'all_pfams', 'total_pfams_in_ecis','total_pfams_NOT_in_eCIS']

df_paths = glob.glob('all_pfam_hits_for_all_genomes/*')
df_total = pd.DataFrame(columns = column_names)

genomes_counted = []
genomes_ecis = []

total_pfams_in_ecis = 0

for path in df_paths:

    df_a = pd.read_table(path)
    
    df_a = df_a.drop(df_a[df_a.gene_oid == 'gene_oid'].index) #clean dataframe of non numerical values
    df_a.gene_oid = df_a.gene_oid.astype('int64')
    df_a['gene_in_eCIS'] = df_a.gene_oid.isin(gene_in_ecis)
    
    df_a = df_a.drop_duplicates(['gene_oid', 'pfam_id'], keep = 'first') #remove multiple pfams from same gene
    
    
    #adding genus and domain to genomes
    df_a = pd.merge(df_a, names, how = 'left', on = 'genome_id') #add Genus for each pfam
    df_a = df_a[df_a.Domain.isin(['Archaea', 'Bacteria', 'Plasmid:Bacteria'])] #take only pfams from these domains
    
    
    # TEST: counting how many pfams found in eCIS operons
    total_pfams_in_ecis += len(df_a[df_a['gene_in_eCIS']])
    print('pfams in eCIS in this df: ', len(df_a[df_a['gene_in_eCIS']]))
    print('total ecis pfam so far: ', total_pfams_in_ecis)
    

    #count how many found in ecis
    pfams_in_ecis = df_a[df_a.gene_in_eCIS == True].pfam_id.value_counts()
    pfams_dict = pfams_in_ecis.to_dict()
    df_a['num_pfams_in_eCIS'] = df_a['pfam_id'].map(pfams_dict)
    df_a.num_pfams_in_eCIS.fillna(0, inplace = True)
    
    #counting all pfam counts
    a = df_a['pfam_id'].value_counts()
    a1 = a.to_dict() #converts to dictionary
    df_a['all_pfams_counts'] = df_a['pfam_id'].map(a1)
    
    
    #counting each pfam occurence not in ecis
    df_a['num_pfam_NOT_in_ecis'] = df_a['all_pfams_counts'] - df_a['num_pfams_in_eCIS']
    
    
    
    df_a['all_pfams'] = len(df_a.pfam_id)
    df_a['total_pfams_in_ecis'] = len(df_a[df_a.gene_in_eCIS == True].pfam_id)
    df_a['total_pfams_NOT_in_eCIS'] = df_a['all_pfams'] - df_a['total_pfams_in_ecis']
    
    
    df_a = df_a.drop_duplicates(subset='pfam_id', keep="last")
    df_a = df_a[['pfam_id', 'pfam_name','num_pfams_in_eCIS', 'all_pfams_counts','num_pfam_NOT_in_ecis', 'all_pfams', 'total_pfams_in_ecis','total_pfams_NOT_in_eCIS']]
    df_total = pd.concat([df_total, df_a])
    
print('finished going throuhg all pfams')



df_total = df_total.astype({'num_pfams_in_eCIS': float,'all_pfams_counts': float,'num_pfam_NOT_in_ecis': float,'all_pfams': float,'total_pfams_in_ecis': float,'total_pfams_NOT_in_eCIS': float})
pfam_names = df_total[['pfam_id', 'pfam_name']].drop_duplicates(subset = 'pfam_id').set_index('pfam_id').to_dict()['pfam_name']


g = df_total.groupby('pfam_id')
df = g.sum()
df = df.reset_index()
df['pfam_name'] = df['pfam_id'].map(pfam_names)
#df = pd.merge(df, pfam_names, how = 'left', on = 'pfam_id')

df['all_pfams'] = df.loc[0].all_pfams
df['total_pfams_in_ecis'] = df.loc[0].total_pfams_in_ecis
df['total_pfams_NOT_in_eCIS'] = df.loc[0].total_pfams_NOT_in_eCIS
 
 
 
df['total_pfams_in_ecis_minus_self'] = df.total_pfams_in_ecis - df.num_pfams_in_eCIS
df['total_pfams_NOT_in_eCIS_minus_self'] = df.total_pfams_NOT_in_eCIS - df.num_pfam_NOT_in_ecis
df.total_pfams_NOT_in_eCIS_minus_self = df.total_pfams_NOT_in_eCIS_minus_self.astype('int64')
 
 
 
 
 # ================================================
 # preparing for fisher exact
 
#preparing matrix for analysis
df['top'] = df[['num_pfams_in_eCIS','num_pfam_NOT_in_ecis']].values.tolist()
df['bottom'] = df[['total_pfams_in_ecis_minus_self', 'total_pfams_NOT_in_eCIS_minus_self']].values.tolist()
 
 
 
index_list = df.index.to_list()

#making dictionary fisher exact for each pfam
FEdict = {}
oddsdict = {}
pvaluedict = {}
 
for i in index_list:
    matrix = np.column_stack((df.loc[i].top, df.loc[i].bottom))
    FE = stats.fisher_exact(matrix)
    FEdict[i] = FE
    oddsdict[i] = FE[0]
    pvaluedict[i] = FE[1]
 
df['odds_ratio'] = df.index.map(oddsdict)
df['pvalue'] = df.index.map(pvaluedict)
 
 
df['FDR_output'] = df['pvalue'].apply(lambda x: multipletests(x, alpha = 0.01, method = 'fdr_bh'))
 
bonf = multipletests(df['pvalue'], alpha = 0.01, method = 'fdr_bh')  
df[['pass_FDR', 'FDR_corrected_pval']] = pd.DataFrame(np.column_stack(bonf[0:2]), index=df.index)
 
 
#add Genus data ofr each pfam
ecis_pfams = pd.read_table('ecis_database_with_pfams_data.tsv') 
ecis_pfams = pd.merge(ecis_pfams, names, how = 'left', on = 'genome_id') #add Genus for each pfam



pfam_genus_in_ecis = ecis_pfams.groupby('pfam_id')['Genus'].nunique().to_dict()


df['number_of_genus_having_this_pfam_in_eCIS_operon'] = df['pfam_id'].map(pfam_genus_in_ecis)
df.number_of_genus_having_this_pfam_in_eCIS_operon.fillna(0, inplace = True)
 
df2 = df[['pfam_id', 'pfam_name', 'number_of_genus_having_this_pfam_in_eCIS_operon', 'odds_ratio',  'FDR_corrected_pval', 'num_pfams_in_eCIS','num_pfam_NOT_in_ecis','all_pfams_counts', 'all_pfams','total_pfams_in_ecis',  'total_pfams_NOT_in_eCIS', 'total_pfams_in_ecis_minus_self', 'total_pfams_NOT_in_eCIS_minus_self','pvalue', 'FDR_output', 'pass_FDR']]
 
df2.to_excel('pfam_enrichment_analysis.xlsx', index = False)


