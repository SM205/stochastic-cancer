#!/bin/bash

#SBATCH --job-name=s4_icfc_permutation     #Name of your job
#SBATCH --cpus-per-task=40    #Number of cores to reserve
#SBATCH --mem-per-cpu=4G     #Amount of RAM/core to reserve
#SBATCH --time=1-00:00:00      #Maximum allocated time
#SBATCH --qos=1day         #Selected queue to allocate your job

ml Python/3.9.6-GCCcore-11.2.0-bare
python ic_fc_s4_permutation.py

