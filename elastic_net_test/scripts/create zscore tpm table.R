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
