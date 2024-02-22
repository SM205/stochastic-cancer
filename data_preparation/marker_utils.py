import copy
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from MarkerCount.marker_count import MarkerCount_Ref, MarkerCount

# Tissue type (311 types) is a narrower classification as compared to tissue class (124 types)

def tissue_type(markers,tissue):
	return markers[markers['tissue_type'] == tissue]

def tissue_class(markers,tissue):
	return markers[markers['tissue_class'] == tissue]

def cell_type(markers,cell):
	return markers[markers['cell_type'] == cell]

def cancer_type(markers,cancer):
	return markers[markers['cancer_type'] == cancer]

def matrix_creation(markers,cell_name="cell_name"):
	markers = markers[[cell_name,'Symbol']]
	markers = markers.dropna()
	markers = markers.groupby(cell_name)['Symbol'].value_counts().unstack('Symbol', fill_value=0).reset_index()
	markers.set_index(cell_name, inplace=True)
	markers[markers > 1] = 1
	return markers

def markercount(adata, marker_matrix):
	adata = adata.transpose()
	X_test = adata.to_df()
	df_res = MarkerCount(X_test, marker_matrix, log_transformed = False, verbose = True )
	adata.obs['cell_type_pred'] = df_res['cell_type_pred'] 
	return adata

def umap(adata,output):
	N_pca = 15
	sc.tl.pca(adata, svd_solver='arpack', n_comps = N_pca)
	sc.pp.neighbors(adata, n_neighbors=10, n_pcs=N_pca)
	sc.tl.umap(adata)
	sc.settings.verbosity = 3
	sc.settings.set_figure_params(figsize=(6, 6), dpi=80, facecolor='white')
	sc.pl.umap(adata, color=['cell_type_pred'], s = 20, legend_fontsize = 12, show = False,save=output)

def classification_results(adata):
	l = []
	for k in adata.obs['cell_type_pred'].drop_duplicates():
		k1 = (adata.obs['cell_type_pred'] == k).sum()
		l.append(k1)
	res = pd.DataFrame(l, index=adata.obs['cell_type_pred'].drop_duplicates())
	return res

def extract_required(adata,to_extract,scRNAseqcounts):
	identity = adata.obs['cell_type_pred']
	final = pd.Series()
	for i in range(len(to_extract)):
		final = final.append(identity[identity == to_extract[i]])
	extraction_index = list(final.index)
	return scRNAseqcounts[extraction_index]