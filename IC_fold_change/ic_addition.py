import pandas as pd
import numpy as np

def get_observed_variance(genes):
    return genes.sum(axis=0).var()

def get_pb_variance(genes):
    fx = genes.sum(axis=1) / genes.shape[1]
    return np.sum([(1-p)*p for p in fx])

def family_filtering_with_ic_calculation(dich):
    ic = []
    filtered_families = []
    for family in family_list:
        genes = dich.reindex(index=familyset[familyset.family_id == family].index)
        mean_gene_per_cell = genes.sum(axis=0).mean()
        non_zero = sum(genes.sum(axis=1) > 0)
        if mean_gene_per_cell > 0.03 and non_zero > 3:
            filtered_families.append(family)
            ic.append(get_observed_variance(genes)/get_pb_variance(genes)) 
    return pd.DataFrame({'family_id': filtered_families, 'IC': ic}).set_index("family_id")

h_dich = pd.read_csv("../Differential_Expression/dichotomised/s2_periphery.csv",index_col="gene_id")
t_dich = pd.read_csv("../Differential_Expression/dichotomised/s2_tumour.csv",index_col="gene_id")
familyset = pd.read_csv("../Differential_Expression/familyset/s2_periphery.csv",usecols=["gene_id", "family_id"],index_col="gene_id")
family_name = pd.read_csv("../clean_pantherhuman.csv")[["family_id","family_name"]].drop_duplicates().set_index("family_id")
family_list = list(familyset["family_id"].drop_duplicates())
h_IC = family_filtering_with_ic_calculation(h_dich)
t_IC = family_filtering_with_ic_calculation(t_dich)
filtered_families = list(t_IC.index.intersection(h_IC.index))

ic_fc = pd.read_csv("ic_fc.csv",index_col="family_id")

ic_fc["IC_tumour"] = t_IC.loc[filtered_families]
ic_fc["IC_healthy"] = h_IC.loc[filtered_families]

ic_fc.to_csv("ic_fc_with_ic.csv")