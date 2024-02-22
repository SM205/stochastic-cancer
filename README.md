# Stochastic Cancer
This set of coddes complements the work "Signatures of cancer cells defined by single-cell gene expression patterns of gene families" (unpublished).  
  
1. [Introduction](#introdction)  
2. [Dependencies](#dependencies) 
3. [Prerequisite](#prerequisites)
4. [Description](#description)  
    - [Data Preparation](#data-preparation)
        - [Cell Selection (MarkerCount) *python script*](#cell-selection)
        - [Salvaging *python script*]()
        - [familyset transformation(Misha) *python script*]()
        - [Dichotomisation(Misha) *python script*]()
    - [Calculation of IC and bootstrapping *python script*]()
    - [Classification of gene families (convert to python script)]()
    - [Contingency table related p-value calculation *R command*]()
    - [Differential gene expression and minP multiple testing correction related *python script*]()
    - [Permutation test based IC fold change calculation *python script*]()
    - [Binomial test and gene ontology enrichment]()
        - [Exclusive families *R notebook*]()
        - [Concurrent families *R notebook*]()
   
# Introdction
  
Different parts of the project were implemented as either standalone python scripts, jupyter notebooks, R scripts or R notebooks. 

Only scripts thought to be imperative for carrying out all the analyses mentioned in the methods section have been provided. 

The following documentation describes all the components in roughly the order in which it occurs in the methods. For the details on the methods and the project itself see the article [pending]().
   
# Dependencies

Python version:
- Python v3.9.12   

Python packages:  
- pandas
- numpy
- 

R version:
- R v
- BioConductor v

R packages:
- tidyverse
- ClusterProfiler
-
  

# Prerequisites
  
clean_pantherhuman.csv - PANTHERv15 Human file provided in the repository. The file was obtained after preprocessing raw panther files 
  
# Description
     
## Data preparation  

Cell selection, familywise transformation and dichotomisation scripts were applied for all the datasets. Salvaging script was only applied for datasets where the mean number of genes detected per cell was less than 2000 (dgpc < 2000).  

### Cell selection  
  


[](https://github.com/d-lowl/stochastic-gene-choice)