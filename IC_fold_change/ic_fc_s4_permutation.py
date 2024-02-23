import pandas as pd
import numpy as np
import statsmodels.stats.multitest
import random
import concurrent.futures
random.seed(218)

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
        if mean_gene_per_cell > 0.1 and non_zero > 3:
            filtered_families.append(family)
            ic.append(get_observed_variance(genes)/get_pb_variance(genes)) 
    return pd.DataFrame({'family_id': filtered_families, 'IC': ic}).set_index("family_id")

def dichotomised_filtering(dich):
    return dich.reindex(index=familyset[familyset.family_id.isin(filtered_families)].index)

def permutation(all_cell_id,number_of_cells_to_select):
    return random.sample(all_cell_id, k=number_of_cells_to_select)

def permutation_ic_fc_calculation(perm_h_columns):
    dich_h = filt_dich[perm_h_columns]
    dich_t = filt_dich.drop(columns=perm_h_columns)
    ic_fc = []
    for family in filtered_families:
        genes_t = dich_t.reindex(index=familyset[familyset.family_id == family].index)
        genes_h = dich_h.reindex(index=familyset[familyset.family_id == family].index)
        ic_t = (get_observed_variance(genes_t)/get_pb_variance(genes_t))
        ic_h = (get_observed_variance(genes_h)/get_pb_variance(genes_h))
        ic_fc.append(ic_t/ic_h)
    return ic_fc

def p_calculation(fc_perm):
    l = []
    for i in range(len(filtered_families)):
        ic_fc = {}
        ic_fc["family_id"] = filtered_families[i]
        ic_fc["ic_fc"] = IC_FC_original.iloc[i].values.item()
        count = sum(abs(ic_fc_perm - fc_perm.iloc[:,i].mean()) >= abs(IC_FC_original.iloc[i] - fc_perm.iloc[:,i].mean()) for ic_fc_perm in fc_perm.iloc[:,i]).values
        pre_p = (count + 1)/(permutations + 1)
        ic_fc["ic_fc_p"] = pre_p.item()
        l.append(ic_fc)
    return pd.DataFrame(l)

permutations = 100000
h_dich = pd.read_csv("s4_periphery/dichotomised_genes.csv",index_col="gene_id")
t_dich = pd.read_csv("s4_tumour/dichotomised_genes.csv",index_col="gene_id")
# familyset of both the tumour/periphery is identical
familyset = pd.read_csv("s4_periphery/familyset.csv",
                        usecols=["gene_id", "family_id"],index_col="gene_id")
family_name = pd.read_csv("clean_pantherhuman.csv")[["family_id","family_name"]].drop_duplicates().set_index("family_id")
family_list = list(familyset["family_id"].drop_duplicates())
h_IC = family_filtering_with_ic_calculation(h_dich)
t_IC = family_filtering_with_ic_calculation(t_dich)
filtered_families = list(t_IC.index.intersection(h_IC.index))
IC_FC_original = (t_IC.loc[filtered_families]/h_IC.loc[filtered_families]).rename(columns={"IC":"IC_FC"})
filt_h_dich = dichotomised_filtering(h_dich).add_suffix("_h")
filt_t_dich = dichotomised_filtering(t_dich).add_suffix("_t")
filt_dich = pd.concat([filt_h_dich,filt_t_dich],axis=1)
permutation_colnames = [permutation(list(filt_dich.columns),filt_h_dich.shape[1]) for _ in range(permutations)]

if __name__ == "__main__":
    with concurrent.futures.ProcessPoolExecutor() as executor:
        fc_perm = executor.map(permutation_ic_fc_calculation,[permutation_colnames[i] for i in range(permutations)])

    fc_perm = pd.DataFrame(fc_perm)
    p_df = p_calculation(fc_perm)
    p_df.set_index("family_id",inplace=True)
    p_df["holm_corrected_ic_fc_p"] = statsmodels.stats.multitest.multipletests(p_df["ic_fc_p"],alpha=0.05,method="holm",is_sorted=False)[1]
    family_name = family_name.loc[filtered_families]
    pd.concat([p_df,family_name],axis=1).to_csv("ic_fc_s4.csv")
