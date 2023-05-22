#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 10:29:05 2023

@author: yaolin
"""



import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.sparse import csr_matrix
import matplotlib as mpl
from matplotlib import pyplot as plt
from ctxcore.genesig import GeneSignature
from pyscenic.aucell import create_rankings, enrichment, aucell, derive_auc_threshold


sc.logging.print_versions()
sc.set_figure_params(frameon=False, figsize=(4, 4), dpi_save=400)

senmayo = pd.read_csv('/Users/yaolin/projects_eriba/Covid_19/SenMayo_human.txt', header=None)
senmayo = senmayo[0].tolist()
print(senmayo)


#############################################   COMBAT  dataset ###################
datafile = "/Users/yaolin/Downloads/COMBAT-CITESeq-DATA.h5ad"
combat = sc.read_h5ad(datafile)
combat

combat.var["feature_types"].value_counts()
combat.obs["Source"].value_counts()
combat.obs["PreExistingHeartDisease"].value_counts()
combat.obs["PreExistingLungDisease"].value_counts()
combat.obs["PreExistingKidneyDisease"].value_counts()
combat.obs["PreExistingDiabetes"].value_counts()
combat.obs["PreExistingHypertension"].value_counts()
combat.obs["PreExistingImmunocompromised"].value_counts()
combat.uns['Source_colors']


ADT = pd.Index.tolist(combat.var_names[20615:])
ADT = [w.replace("AB_", "") for w in ADT]
surface_markers = ['DEP1', 'PTPRJ', 'CD148', 'B2MG', 'B2M', 'CD264', 'TNFRSF10D', 'TRAILR4', 'CD36', 'ICAM-1', 'MDA-VIMENTIN', 'DPP4', 'CD26', 'NOTCH1', 'NOTCH3', 'SCAMP4', 'MICA', 'MICB', 'ULBP2', 'uPAR', 'CD87']
print(set(ADT).intersection(surface_markers))

sc.pl.umap(combat, color='CDKN1A', edges=True)
sc.pl.umap(combat, color='AB_CD36')

combat_covid = combat[(combat.obs['Source'] != 'Flu') & (combat.obs['Source'] != 'Sepsis') & (combat.obs['Source'] != 'COVID_HCW_MILD') & (combat.obs['Source'] != 'COVID_LDN')]
sc.pl.umap(combat_covid, color='CDKN1A')
sc.pl.umap(combat_covid, color='AB_CD36', vmin='p10')
sc.pl.violin(combat_covid, ['CDKN1A', 'AB_CD36'], groupby='Source', rotation=90, stripplot=False )

sc.pl.umap(combat, color='Annotation_major_subset', legend_loc='on data', 
           frameon=False, legend_fontsize=4, legend_fontoutline=1,
           title="Annotation_major_subset", palette='Set1')


sc.pl.umap(combat, color = 'batch', groups = ['batch1'] )
def cluster_small_multiples(adata, clust_key, size=60, frameon=False, legend_loc=None, **kwargs):
    tmp = adata.copy()
    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]
    sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values, color=adata.obs[clust_key].cat.categories.tolist(), size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)
cluster_small_multiples(combat, 'DiseaseClassification')

ax = plt.subplot()
sc.pl.umap(combat, ax=ax, show=False)
sc.pl.umap(combat[combat.obs["DiseaseClassification"].isin(["COVID-19;MONDO:0100096"])], color="AB_CD36", ax=ax, show=False)
plt.show()




######################## monocyte subset #####################
# classical monocytes: expression aprofiling of 'CDKN1A', 'AB_CD36'
combat_cmono = combat[combat.obs['Annotation_major_subset'] == 'cMono']
combat_cmono_covid = combat_cmono[(combat_cmono.obs['Source'] != 'Flu') & (combat_cmono.obs['Source'] != 'Sepsis') & (combat_cmono.obs['Source'] != 'COVID_HCW_MILD') & (combat_cmono.obs['Source'] != 'COVID_LDN')]
sc.pl.violin(combat_cmono_covid, ['CDKN1A', 'AB_CD36'], groupby='Source', rotation=90, stripplot=False )

# non classical monocytes: expression aprofiling of 'CDKN1A', 'AB_CD36'
combat_ncmono = combat[combat.obs['Annotation_major_subset'] == 'ncMono']
sc.pl.violin(combat_ncmono, ['CDKN1A', 'AB_CD36'], groupby='Source', rotation=90)

combat_ncmono_covid = combat_ncmono[(combat_ncmono.obs['Source'] != 'Flu') & (combat_ncmono.obs['Source'] != 'Sepsis') & (combat_ncmono.obs['Source'] != 'COVID_HCW_MILD') & (combat_ncmono.obs['Source'] != 'COVID_LDN')]
sc.pl.heatmap(combat_ncmono, senmayo, groupby='Source', swap_axes=True, cmap='bwr', show_gene_labels=True, dendrogram=True)
sc.pl.violin(combat_ncmono_covid, ['CDKN1A', 'AB_CD36'], groupby='Source', rotation=45, stripplot=False)
sc.pl.violin(combat_ncmono_covid, ['CDKN1A', 'AB_CD26'], groupby='Source', rotation=90)
sc.pl.violin(combat_ncmono_covid, ['AB_CD14', 'AB_CD16'], groupby='Source', rotation=90)




# non classical monocytes: exporting to matrix for analysis by seurat in R 
combat_ncmono_covid.write('combat_ncmono_covid.h5ad')
t=combat_ncmono_covid.X.toarray() # normalized and log transformed
pd.DataFrame(data=t, index=combat_ncmono_covid.obs_names, columns=combat_ncmono_covid.var_names).to_csv('combat_ncmono_covid_normalized.csv')
t=combat_ncmono_covid.layers['raw'].toarray() # raw counts
pd.DataFrame(data=t, index=combat_ncmono_covid.obs_names, columns=combat_ncmono_covid.var_names).to_csv('combat_ncmono_covid_raw.csv')
combat_ncmono_covid.obs.to_csv('combat_obs.csv')
combat_ncmono_covid.var.to_csv('combat_var.csv')
pd.DataFrame(data=combat_ncmono_covid.obsm['X_umap'], index=combat_ncmono_covid.obs_names).to_csv('combat_ncmono_covid_umap.csv')


# senescent genes expression 
combat_ncmono_covid_senmayo = combat_ncmono_covid[:,combat_ncmono_covid.var_names.isin(senmayo)]
sc.pl.heatmap(combat_ncmono_covid_senmayo, combat_ncmono_covid_senmayo.var_names, groupby='Source',standard_scale='var', cmap='bwr', show_gene_labels=True, dendrogram=True)
combat_ncmono_covid_senmayo.copy().write_csvs('combat_ncmono_covid_senmayo.csv', skip_data=False)



sc.set_figure_params(scanpy=True, fontsize=20)
sc.pl.heatmap(combat_ncmono_covid, senmayo, groupby='Source', log=True, standard_scale='var', cmap='bwr', show_gene_labels=True, dendrogram=True)
sc.pl.stacked_violin(combat_ncmono_covid, senmayo, groupby='Source', swap_axes=False, dendrogram=True)
sc.pl.dotplot(combat_ncmono_covid, senmayo, groupby='Source', dendrogram=True)
sc.pl.rank_genes_groups_heatmap(combat_ncmono_covid, groupby='Source', n_genes=10)


help(sc.pp.calculate_qc_metrics)
combat_ncmono_covid
sc.pp.calculate_qc_metrics(combat_ncmono_covid, inplace=True)
sc.pl.violin(combat_ncmono_covid, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True)


