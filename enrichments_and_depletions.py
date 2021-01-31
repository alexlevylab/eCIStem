#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 16:31:56 2020

@author: alexgeller
"""

import pandas as pd
import numpy as np 
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import math
from matplotlib.font_manager import FontProperties


#-------------------
#import ecis database and import all IMG database

pd.set_option('display.max_columns', 55)

db_metadata_all = pd.read_csv(r"~/Desktop/quarantine_work/anchors_expanded_by_10_round2_CDS_only_genomes_only_top50HMMs_filter10HMM_AddedPfamAnnot_filter_eCIS_specific_genes_DD_IMG_METADATA", sep = '\t')
latest_ecis_db = pd.read_csv(r"~/Desktop/quarantine_work/anchors_expanded_by_10_round2_CDS_only_genomes_only_top50HMMs_filter10HMM_AddedPfamAnnot_filter_eCIS_specific_genes_DD_PUBLIC_hmmParams_strand_accessoryAnnot_ScaledBackTo4Accessories_40percentClusters_v2", sep = '\t')

ecis =  latest_ecis_db.merge(db_metadata_all, how = 'left') #add IMG metadata to the final dataframe
ecis.drop_duplicates('genome_ID', inplace = True)


#import all IMG database that are public (in pieces becaused IMG only handles 20,000 at a time in the genome cart)
aa = pd.read_csv(r"~/Desktop/quarantine_work/IMG_aa.csv")
ab = pd.read_csv(r"~/Desktop/quarantine_work/IMG_ab.csv")
ac = pd.read_csv(r"~/Desktop/quarantine_work/IMG_ac.csv")
ad = pd.read_csv(r"~/Desktop/quarantine_work/IMG_ad.csv")
ae = pd.read_csv(r"~/Desktop/quarantine_work/IMG_ae.csv")
af = pd.read_csv(r"~/Desktop/quarantine_work/IMG_af.csv")


db = pd.concat([aa, ab, ac, ad, ae, af], ignore_index = True, names = ['taxon_oid', 'Domain',
       'Genome Name / Sample Name',
       'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain',
       'Biotic Relationships', 'Ecosystem', 'Ecosystem Category',
       'Ecosystem Subtype', 'Ecosystem Type', 'Ecotype', 'Habitat',
       'Host Name', 'Isolation', 'Sample Body Site', 'Sample Body Subsite',
       'Specific Ecosystem'])


db = db.loc[db['Domain'].isin(['Bacteria', 'Archaea', 'Plasmid:Bacteria'])]

#------------------


# function to fill in parts of table that will be used for Fisher exact test
def enr(df):
    truetrue = df[taxon + '_eCIS']
    truefalse = df[taxon + '_database'] - df[taxon + '_eCIS']
    return [truetrue, truefalse]


#------------------


taxa = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Ecosystem', 'Ecosystem Category', 'Ecosystem Subtype', 'Ecosystem Type', 'Sample Body Site', 'Sample Body Subsite',  'Habitat',
       'Host Name', 'Isolation', 'Biotic Relationships', 'Motility', 'Oxygen Requirement', 'Phenotype', 'Temperature Range']

#----


outputdf = pd.DataFrame()


#----    
for taxon in taxa: # loop through each category under question, e.g. Phylum. Can do taxa[0:1] to test one at a time. Some are not taxa, but are called taxa anyway
    #-----------------
    # Reset outdf to be a blank dataframe for new loop instance
    
    out_df = pd.DataFrame({'taxon': [], # Reset outdf to be a blank dataframe
                    'group': [],
                    'odds_ratio': [],
                    'pvalue': [],
                    'significant': []})    

    #-----------------
    
    # calculate number of appearances (frequency) of each group, e.g. Proteobacteria in ecis database vs IMG database
    ecis_freq = ecis[taxon].value_counts().reset_index()
    database_freq = db[taxon].value_counts().reset_index()
    
    #merge frequences
    merged = pd.merge(database_freq, ecis_freq, how = 'left', on = 'index', suffixes=('_database', '_eCIS'))
    
    #clean up table
    merged['percent'] = merged[taxon + '_eCIS'] / merged[taxon + '_database'] 
    merged.sort_values([(taxon + '_database'), 'percent'])
    merged['Taxon'] = merged['index']
    merged2 = merged[['Taxon', taxon + '_database', taxon + '_eCIS']]
    
    #remove unclassified genomes
    merged2= merged2[~(merged2['Taxon'].isin(['unclassified', 'Unclassified', 'Placeholder']))]
    
    #final product is merged2
    merged2.fillna(value = 0.0, inplace = True)
    
    print(merged2) 
    
    #-----------------
    # go through each group (e.g. Proteobacteria, Actinobacteria etc.)
    
    for group in merged2.Taxon.unique():
        print("GROUP= ", group)
        
        #fill in table for Fisher Exact Test
        top = merged2[merged2['Taxon'] == group].apply(enr, axis = 1)
        part2 = merged2[merged2['Taxon'] != group]
        bottom_truetrue = part2[taxon + '_eCIS'].sum()
        bottom_truefalse = (part2[taxon + '_database'] - part2[taxon + '_eCIS']).sum()
        bottom = [bottom_truetrue, bottom_truefalse]
        matrix = np.column_stack((top.values[0], bottom))
        print(matrix)
        
        # calculate fisher exact test for each group; e.g. are eCIS genomes enriched in proteobacteria?
        oddsratio, pvalue = stats.fisher_exact(matrix)
        print(oddsratio)
        print(pvalue)
        
        #add results to out_df one by one, so at the end all groups are together in one table, e.g proteobacteria fisher test, actinobacteria fisher test... 
        if pvalue < 0.01:
            print("-----SIGNIFICANT^^^---------")
            out_df = out_df.append(pd.Series([taxon, group, oddsratio, pvalue, True], index=['taxon','group', 'odds_ratio', 'pvalue', 'significant']), ignore_index=True)
        else:
            print("not significant :( :( :( :( ")
            out_df = out_df.append(pd.Series([taxon, group, oddsratio, pvalue, False], index=['taxon','group', 'odds_ratio', 'pvalue', 'significant']), ignore_index=True)



    merged2['percent'] =  (merged2[taxon + '_eCIS'] / merged2[taxon + '_database']) *100

#    ----------------------
    
    final = pd.merge(merged2, out_df, how = 'left', left_on = 'Taxon', right_on = 'group')
    
    final.fillna(value = False, inplace = True)
#    print(final)
    
    final2 = final
    
    print(final)
    
    
#---------------------------------------------------------------------------------------------------
#     # To export final product: 
#     
    outtocsv = pd.DataFrame()
    outtocsv = final
    outtocsv['num_database'] = final[taxon + '_database']
    outtocsv['num_ecis'] = final[taxon + '_eCIS']
    outtocsv['category'] = taxon
    outtocsv = outtocsv[['category', 'Taxon', 'num_database', 'num_ecis', 'percent', 'odds_ratio', 'pvalue']]

    outputdf = pd.concat([outputdf, outtocsv])
    outtocsv = pd.DataFrame()
    


bonf = multipletests(outputdf['pvalue'], alpha = 0.01, method = 'fdr_bh')  
outputdf[['pass_FDR', 'FDR_corrected_pval']] = pd.DataFrame(np.column_stack(bonf[0:2]), index=outputdf.index)


print(outputdf)

outputdf.to_csv(r"~/Desktop/quarantine_work/new_database_enrichments_FIXED_FDR_fixedDBnumber_GENUSONLY.csv", index = False)
outputdf.loc[outputdf['pass_FDR']== True].to_csv(r"~/Desktop/quarantine_work/new_database_enrichments_statistically_significant_fixedDBnumber_GENUSONLY.csv", index = False)    
    #----------------

 
    
    
    
    
    
    
    
    

