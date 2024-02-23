import glob
import pandas as pd

path = "study/final_data/*.csv"

for file in glob.glob(path):
	df = pd.read_csv(file)
	df = df.rename(columns={'Unnamed: 0':'gene_id'})
	df = df.set_index("gene_id")
	df.to_csv(file)