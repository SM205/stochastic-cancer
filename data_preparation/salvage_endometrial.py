import pandas as pd

def salvage(df):
    while (df>0).sum().quantile(0.1) < 2000:
        del df[(df>0).sum().sort_values().index[0]]
    return df

e3 = pd.read_csv("filtered_EC3.csv",index_col="gene_id")
salvage(e3).to_csv("salvaged_EC3.csv")