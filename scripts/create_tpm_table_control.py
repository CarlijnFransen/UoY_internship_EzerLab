#!/bin/python

#script to convert RSEM output in various folders to one TPM table
from os import listdir
import sys
import pandas as pd

path = "rsem_output/"
sample_folders = sorted(listdir(path))
#print(sample_folders)
complete_tpm_table_bool = True
sample_names_list = [""]
for sample in sample_folders:
    if sample.startswith("."):
        pass
    elif sample.endswith("T"):
        pass   
    elif sample.endswith("_T_concatenated_4"):
        pass
    else:
        if sample.strip().endswith("_concatenated_4"):
            short_sample = sample.strip("_concatenated_4")
        else: 
            short_sample = sample
        short_sample_no_c = short_sample.strip("_C")
        sample_names_list.append(short_sample_no_c)
print(sorted(sample_names_list))
sorted_sample_names_list = sorted(sample_names_list)
for short_sample_no_c in sorted_sample_names_list:
    if short_sample_no_c == "":
        continue
    short_sample = short_sample_no_c+"_C"
    sample = short_sample
    if short_sample_no_c ==  "WT_63":
        continue
    elif short_sample_no_c == "WT_35":
        continue
    try:
        TPM_path = path + sample + "/"+ short_sample+".genes.results"
        TPM_file = pd.read_csv(TPM_path , sep='\t')
    except:
        sample = short_sample_no_c+"_C_concatenated_4"
        TPM_path = path + sample + "/"+ short_sample+".genes.results"
        TPM_file = pd.read_csv(TPM_path , sep='\t')

    if complete_tpm_table_bool == True:
        complete_tpm_table_bool = False
        complete_tpm_table = TPM_file.filter(["gene_id","TPM"])
        complete_tpm_table.rename(columns = {"gene_id":""}, inplace=True)
        complete_tpm_table.rename(columns = {"TPM":short_sample_no_c}, inplace=True)
    else:
        complete_tpm_table[short_sample_no_c] = TPM_file["TPM"]
        
      
complete_tpm_table.to_csv('tnz_ws2_ril_TPM_table_only_control_2.csv',  index=False, header=True)
