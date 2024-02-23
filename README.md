# Stochastic Cancer
This set of coddes complements the work "Signatures of cancer cells defined by single-cell gene expression patterns of gene families" (unpublished).  
  
1. [Introduction](#introduction)  
2. [Dependencies](#dependencies) 
3. [Prerequisite](#prerequisites)
4. [Description](#description)  
    - [Data Preparation](#data-preparation)
        - [Cell Selection (MarkerCount)](#cell-selection)
        - [Cell filtering *R script*]()
        - [Salvaging *python script*](#salvaging)
        - [familyset transformation](#familyset-transformation)
        - [Dichotomisation](#dichotomisation)
    - [Calculation of IC and bootstrapping](#calculation-of-ic-and-bootstrapping)
    - [Classification of gene families](#gene-family-classification)
        - [Exclusive](#exclusive)
        - [Concurrent](#concurrent)
    - [Differential gene expression and minP multiple testing correction](#differential-gene-expression-and-minp-multiple-testing-correction)
        - [TMM normalisation](#tmm-normalisation)
        - [Wilcoxon test and minP adjustment](#wilcoxon-test-and-minp-adjustment)
    - [Permutation test based IC fold change calculation](#permutation-test-based-ic-fold-change-calculation)
    - [Binomial test and gene ontology enrichment](#binomial-test-and-go-enrichment)
        - [Exclusive families](#for-exclusive-families)
        - [Concurrent families](#for-concurrent-families)
   
# Introduction
  
Different parts of the project were implemented as either standalone python scripts, jupyter notebooks, R scripts or R notebooks. 

Only scripts thought to be imperative for carrying out the analyses mentioned in the methods section have been provided. (scripts used to generate some of the figures have not been provided, but can be provided if required)

The following documentation describes all the components in roughly the order in which it occurs in the methods. For the details on the methods and the project itself see the article [pending]().
   
# Dependencies

Python version:
- Python v3.9.12   

Python packages:
- os  
- glob
- pandas
- numpy
- MarkerCount (markercount)
- multiprocessing (minP)
- random (minP,icfc)
- scipy (minP)
- statsmodels (icfc)
- concurrent.futures (icfc)

R version:
- R v4.3.2
- BioConductor v3.17

R packages:
- tidyverse
- edgeR (tmm)
- readxl
- ClusterProfiler (GO)
- org.Hs.eg.db (GO)
  

# Prerequisites
  
clean_pantherhuman.csv - PANTHERv15 Human file provided in the repository. The file was obtained after preprocessing raw panther files. This was obtained using the script preprocess_panther.py from the [repository from our previous publication](https://github.com/d-lowl/stochastic-gene-choice).
  
Cell_marker_Human.xlsx - This is a list of all marker genes in humans from [CellMarker 2.0](http://bio-bigdata.hrbmu.edu.cn/CellMarker/), and a different subset has been used for the different types of cancers for classifying the cells with [MarkerCount](https://pubmed.ncbi.nlm.nih.gov/35782735/)  
  
# Description
     
## Data preparation  

Cell selection, cell filtering, familywise transformation and dichotomisation scripts were applied for all the datasets. Salvaging script was only applied for datasets where the mean number of genes detected per cell was less than 2000 (dgpc < 2000).  

There is a manual step in both cell selection and cell filtering which means that the script has minor differences between each study from where the data was extracted. The script uploaded have had dataset/study specific information removed. Some require manual changes (like marker information, cell type selection information, etc.) to work.

- ### **Cell selection**
  
    Two scripts are involved in this step, ***marker_count.py*** and ***marker_utils.py***. The marker_utils script contains various functions used by marker_count script. marker_utils also imports MarkerCount and has functions for all the analysis. The marker_count script provides the criteria used for obtaining the marker genes to be used in further analysis (list provided in Table S1) and outputs the file necessary (in a Predictions folder that needs to be created) for deciding upon the manual filtering ([sample]_count.csv) and the file necessary for manual filtering ([sample]_table.csv) using the cell filtering R script.

    Input - [sample].csv, Cell_marker_Human.xlsx, information on markers to be used need to be manually defined in the script along with the path where all the samples from a study can be found (provided in Table S1)

    Output - Predictions/[sample]_count.csv and Predictions/[sample]_table.csv 

  
- ### **Cell filtering**

    The ***cell_filtering.R*** script used for cell filtering is different for each dataset, because the list of cell types to retain varies between different studies. Due to this, the set has to be manually inserted to the script to obtain the filtered versions of the dataset. Besides that, the path to the study folder also needs to be added manually, after which it can work for all samples in a study as long as enough variables are created in the same pattern. (another python script to add gene_id)

    Input - [sample].csv, Predictions/[sample]_count.csv, Predictions/[sample]_table.csv 
and information about cell types to select  
      
    Output -  filtered_sample.csv (one for each sample)
      
    To make the output in a form that is readable by transform_familyset, the python script ***adding_gene_id.py*** was used.


- ### **Salvaging**
    ***Salvage.py*** removes cells with the lowest gene counts until mean gene per cell becomes > 2000. Only applied to datasets with > 100 cells but mean gene per cell < 2000 after filtering. There is a need to check the cell number after salvaging to know if the sample can be retained for further anaylsis or not.

    Input - filtered_sample.csv
    Output - salvaged_sample.csv
  
- ### **familyset transformation**
    ***transform_familyset.py*** script adapted from [previous publication](https://github.com/d-lowl/stochastic-gene-choice) with minor changes to obtain a familyset from each dataset.
  
    Input - clean_pantherhuman.csv, [filtered/salvaged]_[sample].csv
      
    Output - familyset.csv

- ### **Dichotomisation**
    ***dichotomise_families.py*** script adapted from [previous publication](https://github.com/d-lowl/stochastic-gene-choice) with no change to obtain the dichotomised gene file as well as a file for familywise thresholds.
      
    Input - familyset.csv, [sample]_[filtered/salvaged].csv, "GTME"/"FM" dichotomisation method
      
    Output - dichotomised_genes.csv, family_thresholds.csv

## Calculation of IC and bootstrapping 
  
The ***features_and_bootstrap.py*** script takes in the dichotomised expression of genes and outputs the bootstrapped IC values along with other information related to the family, including mean gene per cell and non zero genes, both of which are important for selection of families, as mentioned in the methods section.

Input - dichotomised_genes.csv, clean_pantherhuman.csv, familyset.csv (familyset and dichotomised_genes of the sample the script is run for)
    
Output - Final_IC_bootstrap.xlsx

  
## Gene family classification  
    Directory structure:  
    * Dataset_folder
        * clean_pantherhuman.csv  
        * Study1
            * Bootstrap  
                * Patient1
                    * features_and_bootstrap.py
                    * Final_IC_bootstrap.xlsx
                    * dichotomised_genes.csv
                    * familyset.csv
                * Patient2
                * ...
        * Study2
        * ...
- ### **Exclusive**
    ***Exclusive.ipynb*** is the python notebook to create an overview of all the families found to be exclusive in at least one sample.  

    Input-  clean_pantherhuman.csv, list of Study1,Study2,... 
      
    Output-  overview_exclusive_table.xlsx (will be further annotated before use in [overrepresentation/binomial test](#for-exclusive-families))

- ### **Concurrent**
    ***Concurrent.Rmd*** is the R notebook to create an overview of all the families found to be concurrent in at least one sample. The directory structure is copied from *Exclusive.ipynb* output so that the patient order is consistent between the two.
  
    Input-  clean_pantherhuman.csv, overview_exclusive_table.xlsx
      
    Output-  concurrent_overview_table.csv (will be further annotated before use in [overrepresentation/binomial test](#for-concurrent-families))  
          
## Differential gene expression and minP multiple testing correction
  
- ### **TMM normalisation**
    TPM values are converted to TMM normalised TPM values using the R package edgeR in the script ***TMM_normalisation.Rmd***. Since there are a total of four files that need to be converted, all of them have been normalised in the same R notebook.     
  
    Input - s2_periphery.csv, s2_tumour.csv, s4_periphery.csv, s4_tumour.csv (All are TPM data files for the glibolastoma dataset)  
  
    Output - tmm_s2_periphery.csv, tmm_s2_tumour.csv, tmm_s4_periphery.csv, tmm_s4_tumour.csv
  
- ### **Wilcoxon test and minP adjustment**
    ***minP_s2.py*** and ***minP_s4.py*** are the two scripts for performing the Wilcoxon test and minP adjustment on files from S2 and S4 respectively. Besides the difference in the input files, the two scripts are identical. 10,000 permutations are performed for minP.

    Input (all from output of TMM normalisation):  
    *minP_s2.py* - tmm_s2_periphery.csv, tmm_s2_tumour.csv  
    *minP_s4.py* - tmm_s4_periphery.csv, tmm_s4_tumour.csv

    Output:   
    *minP_s2.py* - s2_DE_permutation_wilcoxon.csv  
    *minP_s4.py* - s4_DE_permutation_wilcoxon.csv

## Permutation test based IC fold change calculation
  
***ic_fc_s2_permutation.py*** and ***ic_fc_s4_permutation.py*** are the two scripts for performing the permutation on the healthy vs cancer cells of S2 and S4 respectively. Besides the difference in the input files, the two scripts are identical. 100,000 permutations are perfomed and bonferroni-holm multiple testing correction is applied.

Directory structure of the patient folder for code execution:  
* patient  
    * ic_fc_patient_permutation.py
    * clean_pantherhuman.csv
    * patient_periphery  
        * familyset.csv
        * dichotomised_genes.csv
    * patient_tumour 
        * familyset.csv
        * dichotomised_genes.csv

Input (The outputs are from *transform_familyset.py* and *dichotomise_families.py*):  
*ic_fc_s2_permutation.py* - s2_periphery/dichotomised_genes.csv, s2_tumour/dichotomised_genes.csv, s2_periphery/familyset.csv, clean_pantherhuman.csv  
*ic_fc_s4_permutation.py* - s4_periphery/dichotomised_genes.csv, s4_tumour/dichotomised_genes.csv, s4_periphery/familyset.csv, clean_pantherhuman.csv

Output:   
*ic_fc_s2_permutation.py* - ic_fc_s2.csv  
*ic_fc_s4_permutation.py* - ic_fc_s4.csv

  
## Binomial test and GO enrichment
  
- ### **For exclusive families**
    ***Exclusive_Binomial_test_analysis.Rmd*** is the R notebook used for binomial testing of exclusive families for overrepresentation and GO enrichment analysis.  
    
    Input - overview_exclusive_table.xlsx (obtained after manual annotation of output from [gene family classification](#gene-family-classification), provided in the binomial test folder), clean_pantherhuman.csv 
      
    Output - As necessary
  
- ### **For concurrent families**
    ***Concurrent_Binomial_test_analysis.Rmd*** is the R notebook used for binomial testing of concurrent families for overrepresentation and GO enrichment analysis.  
    
    Input - overview_concurrent_table.xlsx (obtained after manual annotation of output from [gene family classification](#gene-family-classification), provided in the binomial test folder), clean_pantherhuman.csv
      
    Output - As necessary