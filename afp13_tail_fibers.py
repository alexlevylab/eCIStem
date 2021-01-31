import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
pd.set_option('display.max_columns', 25)


ad_input = 'afp13_GENOME_GENE_TOP75PERCENTILE_INCLUSIVELIST629members_out6_10239_minus_minus_again'

ph_input = 'afp13_GENOME_GENE_TOP75PERCENTILE_INCLUSIVELIST629members_out6_28883'



'''
this function plots blast hits for eCIS fiber tails vs. either Viruses or Tailed phages and compares hits

outfmt is: 6 + (extras listed here)
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle sstart send
'''

ad = pd.read_csv(r"/Users/alexgeller/Desktop/quarantine_work/" + ad_input, sep = "\t", header = None,
                 names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"])

ph = pd.read_csv(r"/Users/alexgeller/Desktop/quarantine_work/" + ph_input, sep = "\t", header = None,
                  names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"])


ad = ad.loc[ad['evalue'] < 0.001]
ph = ph.loc[ph['evalue'] < 0.001]


names_ad = pd.read_csv("/Users/alexgeller/Desktop/quarantine_work/virus_names_aug2020", header = None,
                       names = ['sseqid', 'subject_organism', 'subject_protein_name'])


names_ph = pd.read_csv("/Users/alexgeller/Desktop/quarantine_work/phage_names_aug2020", header = None,
                       names = ['sseqid', 'subject_organism', 'subject_protein_name'])


ad = ad.merge(names_ad, on = 'sseqid', how = 'left')
ph = ph.merge(names_ph, on = 'sseqid', how = 'left')


ad_max = ad.sort_values(['qseqid', 'bitscore'], ascending = False).drop_duplicates('qseqid', keep = 'first')
ph_max = ph.sort_values(['qseqid', 'bitscore'], ascending = False).drop_duplicates('qseqid', keep = 'first')

print(ad_max)
print(ph_max)


merged = pd.merge(ad_max, ph_max, on = 'qseqid', suffixes = ('_virus', '_phage'), how = 'outer')

merged.fillna(0, inplace = True)

merged['difference_in_score'] =  merged['bitscore_virus'] -  merged['bitscore_phage'] 
merged['abs'] = abs(merged['difference_in_score'])

merged.sort_values("difference_in_score", inplace = True)

#---------
fig, ax = plt.subplots(figsize = (25, 10))

ax.bar(x = merged['qseqid'].astype(str), height = merged['difference_in_score'])
plt.xticks(fontsize=5, rotation=90)
ax.set_ylabel("Top Virus Bitscore - Top Phage Bitscore")
ax.set_xlabel("afp13 gene ID")

#plt.show()
#plt.savefig("afp13_graph_best_viral_hit_bitscore_minus_best_phage_hit_bitscore.png")
#---------



def score(df):
   if df['abs'] < 20:
       return "no clear winner"
   
   elif df['difference_in_score'] < 0:
       return "phage"

   elif df['difference_in_score'] > 0:
       return "virus"
   
   else:
       return "ERROR"


merged['winner'] = merged.apply(score, axis = 1)

#merged.to_csv('/Users/alexgeller/Desktop/quarantine_work/afp13_winners_list', index = False)



