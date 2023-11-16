#!/bin/python

#script to convert RSEM output in various folders to one TPM table
from os import listdir
import sys
import pandas as pd

path = "rsem_output/"
sample_folders = listdir(path)
complete_tpm_table_bool = True
sample_names_list = ["gene_id"]
for sample in sample_folders:
    if sample.startswith("."):
        continue
    else:
        if sample.strip().endswith("_concatenated_2"):
            short_sample = sample.strip("_concatenated_2")
        else: 
            short_sample = sample
        sample_names_list.append(short_sample)
        #print(short_sample)
        TPM_path = path + sample + "/"+ short_sample+".genes.results"
        TPM_file = pd.read_csv(TPM_path , sep='\t')
        if complete_tpm_table_bool == True:
            complete_tpm_table_bool = False
            complete_tpm_table = TPM_file.filter(["gene_id","TPM"])
            complete_tpm_table.rename(columns = {"TPM":short_sample}, inplace=True)
        else:
            complete_tpm_table[short_sample] = TPM_file["TPM"]
        
        
complete_tpm_table.to_csv('tnz_ws2_ril_TPM_table.csv',  index=False, header=True)