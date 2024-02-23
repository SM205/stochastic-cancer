import pandas as pd

def salvage(df):
    while (df>0).sum().quantile(0.1) < 2000:
        del df[(df>0).sum().sort_values().index[0]]
    return df

sample = pd.read_csv("filtered_sample.csv",index_col="gene_id")
salvage(sample).to_csv("salvaged_sample.csv")