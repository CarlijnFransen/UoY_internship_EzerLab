TPM_table_only_control <- read.csv("~/WUR/Internship/Internship/to_git/UoY_internship_EzerLab/elastic_net_test/data/tnz_ws2_ril_TPM_table_only_control_2.csv", row.names=1)


zscore <- function(mat, direction=1){
  apply(mat,direction,  function(i){
    (i-mean(i))/sd(i)
  })
}
library(dplyr)



zscore_tpm_table <- zscore(TPM_table_only_control)
new_zscore_tpm_table <- zscore_tpm_table[ , colSums(is.na(zscore_tpm_table))==0]
write.csv(new_zscore_tpm_table, "~/WUR/Internship/Internship/to_git/UoY_internship_EzerLab/elastic_net_test/data/zscore_tpm_table_only_control.csv", row.names = T, sep = "\t")

zscore_transposed <- as.data.frame(t(new_zscore_tpm_table))
########################################################################
########################### Filter stuff ###############################
########################################################################

# filter for transcription factors
Ath_TF_list <- read.delim("~/WUR/Internship/Internship/to_git/UoY_internship_EzerLab/elastic_net/data/Ath_TF_list.txt", row.names = 1)

TF_filtered_list <- unique(Ath_TF_list$Gene_ID)

TPM_table_only_control_filtered_tf <- TPM_table_only_control[rownames(TPM_table_only_control) %in% TF_filtered_list,]
zscore_tpm_table_filtered_tf <- zscore_transposed[rownames(zscore_transposed) %in% TF_filtered_list,]

write.csv(TPM_table_only_control_filtered_tf, "~/WUR/Internship/Internship/to_git/UoY_internship_EzerLab/elastic_net_test/data/tpm_table_filterd_tf.csv", row.names = T)
write.csv(zscore_tpm_table_filtered_tf, "~/WUR/Internship/Internship/to_git/UoY_internship_EzerLab/elastic_net_test/data/zscore_tpm_table_filterd_tf.csv", row.names = T)



#export this to a new file
TPM_table_only_control_filtered <- TPM_table_only_control[rownames(TPM_table_only_control) %in% unique_combined_selected_samples,]
write.csv(TPM_table_only_control_filtered, "~/WUR/Internship/Internship/to_git/UoY_internship_EzerLab/elastic_net_test/data/tpm_table_filterd_quadrants.csv", row.names = T)

# also filter z_scores_tpm table
#zscore_transposed <- as.data.frame(t(zscore_tpm_table))
zscore_tpm_table_filtered <- zscore_transposed[rownames(zscore_transposed) %in% unique_combined_selected_samples,]
write.csv(zscore_tpm_table_filtered, "~/WUR/Internship/Internship/to_git/UoY_internship_EzerLab/elastic_net_test/data/zscore_tpm_table_filterd_quadrants.csv", row.names = T)
 
# get sample color
condition_rowname_2 <- rownames(condition_subdf[rownames(condition_subdf) %in% rownames(condition_df[["q2"]]),]) #1
condition_rowname_3 <- rownames(condition_subdf[rownames(condition_subdf) %in% rownames(condition_df[["q3"]]),]) #2
samples_rowname_2 <- rownames(samples_subdf[rownames(samples_subdf) %in% rownames(samples_df[["q2"]]),]) #3
samples_rowname_3 <- rownames(samples_subdf[rownames(samples_subdf) %in% rownames(samples_df[["q3"]]),]) #4

TPM_table_rowname <- rownames(TPM_table_only_control)

color_vector <- vector("string", length=28775)

# none        | 0
# condition 2 | 1
# condition 3 | 2
# sample    2 | 3
# sample    3 | 4
# C2 + S2     | 5
# C2 + S3     | 6
# C3 + S2     | 7
# C3 + S3     | 8

combi_vect <- c()

for (i in TPM_table_rowname){
  if (i %in% condition_rowname_2){
    if (i %in% samples_rowname_2){
      combi_vect <- append(combi_vect, 5)
    }
    else if (i %in% samples_rowname_3){
      combi_vect <- append(combi_vect, 6)
    }
    else{
      combi_vect <- append(combi_vect, 1)
    }
  }
  else if (i %in% condition_rowname_3){
    if (i %in% samples_rowname_2){
      combi_vect <- append(combi_vect, 7)
    }
    else if (i %in% samples_rowname_3){
      combi_vect <- append(combi_vect, 8)
    }
    else{
      combi_vect <- append(combi_vect, 2)
    }
  }
  else{
    if (i %in% samples_rowname_2){
      combi_vect <- append(combi_vect, 3)
    }
    else if (i %in% samples_rowname_3){
      combi_vect <- append(combi_vect, 4)
    }
    else{
      combi_vect <- append(combi_vect, 0)
    }
  }
  
}

color_df <- data.frame(TPM_table_rowname, combi_vect , row.names = 1 )

write.csv(color_df, "~/WUR/Internship/Internship/to_git/UoY_internship_EzerLab/elastic_net/data/color_tpm.csv", row.names = T)


