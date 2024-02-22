import glob
import pandas as pd
import scanpy as sc
import marker_utils

path = "/Users/mondal0000/Documents/CellMarkerDB_V2/GSE99330/*.csv"

cell_markers = pd.read_excel('/Users/mondal0000/Documents/CellMarkerDB_V2/Marker_Data/Cell_marker_Human.xlsx')
markers = marker_utils.cell_type(cell_markers,"Cancer cell")
markers = marker_utils.cancer_type(markers,"Melanoma")
marker_matrix = marker_utils.matrix_creation(markers, "cell_name")
marker_matrix = marker_matrix[marker_matrix.sum(axis=1) > 4]

# For some reason, there is a need to create the directory Predictions and it cannot be created by python on its own

for file in glob.glob(path):
	split_path = file.split("/")
	file_name = split_path[-1][:-4]
	adata = sc.read_csv(file)
	adata = marker_utils.markercount(adata,marker_matrix)
	marker_utils.classification_results(adata)
	marker_utils.umap(adata,f'_{file_name}.png')
	l = []
	for k in adata.obs['cell_type_pred'].drop_duplicates():
		k1 = (adata.obs['cell_type_pred'] == k).sum()
		l.append(k1)
	res = pd.DataFrame(l, index=adata.obs['cell_type_pred'].drop_duplicates())
	res.sort_values(by=0, ascending=False).to_csv(f'Predictions/{file_name}_count.csv')
	pd.DataFrame(adata.obs).to_csv(f'Predictions/{file_name}_table.csv')
