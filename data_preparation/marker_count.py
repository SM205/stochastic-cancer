import glob
import pandas as pd
import scanpy as sc
import marker_utils

path = study_accession/*.csv #path to folder containing sample files in csv format

cell_markers = pd.read_excel('Cell_marker_Human.xlsx')
markers = marker_utils.cell_type(cell_markers,Cell Type)
# usually only one of the two types below will be used
# In certain studies with multiple criteria, another step is added
# Example given below
# markers_1 = markers[markers['cancer_type'].str.contains("Triple-Negative Breast Cancer")]
# markers_2 = markers[markers['cancer_type'].str.contains("Triple-negative breast cancer")]
# markers = pd.concat([markers_1,markers_2])
markers = marker_utils.tissue_type(markers,Tissue Type) 
markers = marker_utils.cancer_type(markers,Cancer Type)

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
