#install.packages("pheatmap")
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install(version = "3.18")


library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(BiocManager)

#Read files for each of the parents and conditions
TNZ_1_C_genes <- read_delim("new_tpm_tables/TNZ_1_C.genes.results", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
TNZ_1_T_genes <- read_delim("new_tpm_tables/TNZ_1_T.genes.results", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
WS_2_C_genes <- read_delim("new_tpm_tables/WS_2_C.genes.results", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
WS_2_T_genes <- read_delim("new_tpm_tables/WS_2_T.genes.results", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

#create one big tpm table and change the column names
tpm_table <- TNZ_1_C_genes[,1]
tpm_table <- cbind(tpm_table, TNZ_1_C_genes[,6])
tpm_table <- cbind(tpm_table, TNZ_1_T_genes[,6])
tpm_table <- cbind(tpm_table, WS_2_C_genes[,6])
tpm_table <- cbind(tpm_table, WS_2_T_genes[,6])

colnames(tpm_table) <- c("genes" ,"TNZ_1_C", "TNZ_1_T","WS_2_C","WS_2_T")

#transpose tpm table
t_tpm_table <- t(tpm_table)
colnames(t_tpm_table) <- t_tpm_table[1,]
t_tpm_table <- t_tpm_table[-1, ] 

rownames(tpm_table) <- tpm_table[,1]
tpm_table[,1] <- NULL


# filter on minimum value of X for all variables and at least one value of 20 in each row
# we discussed all numbers minimum value of 5
filtered_tpm_table <- tpm_table %>% filter_all(all_vars(.>=5)) %>% filter_all(any_vars(.>=20))

##### PLOTS #####
#change number of plots in figure to 2X2
par(mfrow=c(2,2))

#plot conditions and ecotypes againts each other again
plot(log(filtered_tpm_table$TNZ_1_C), log(filtered_tpm_table$TNZ_1_T),xlab="TNZ control",ylab = "TNZ light"
     #, main= "Comparison of log TPM between TNZ control and TNZ light conditions"
     , pch=20, col=rgb(0.1,0.1,0.1,0.3))
abline(c(0,1), col = "blue", lwd=2, untf=TRUE)

plot(log(filtered_tpm_table$WS_2_C), log(filtered_tpm_table$WS_2_T),xlab="WS-2 control",ylab = "WS-2 light"
     #, main= "Comparison of log TPM between WS-2 control and WS-2 light conditions"
     , pch=20, col=rgb(0.1,0.1,0.1,0.3))
abline(c(0,1), col = "blue", lwd=2, untf=TRUE)

plot(log(filtered_tpm_table$WS_2_C), log(filtered_tpm_table$TNZ_1_C),xlab="WS-2 control",ylab = "TNZ control"
     #, main= "Comparison of log TPM between WS-2 control and TNZ control conditions"
     , pch=20, col=rgb(0.1,0.1,0.1,0.3))
abline(c(0,1), col = "blue", lwd=2, untf=TRUE)

plot(log(filtered_tpm_table$WS_2_T), log(filtered_tpm_table$TNZ_1_T),xlab="WS-2 light",ylab = "TNZ light"
     #, main= "Comparison of log TPM between WS-2 light and TNZ light conditions"
     , pch=20, col=rgb(0.1,0.1,0.1,0.3))
abline(c(0,1), col = "blue", lwd=2, untf=TRUE)


par(mfrow=c(1,1)) #change again to one graph per figure
# plot the scatterplot for ecotypes
plot(log(filtered_tpm_table$TNZ_1_C/filtered_tpm_table$WS_2_C)
     , log(filtered_tpm_table$TNZ_1_T/filtered_tpm_table$WS_2_T) , xlab="log fold change TNZ and WS-2 control "
     , ylab = "log fold change TNZ and WS-2 light"
     #, main= "Comparison of log TPM between WS-2 light and TNZ light conditions"     
     , pch=20, col=rgb(0.1,0.1,0.1,0.3))
abline(v=0,col='blue', lwd=2)
abline(h=0,col='blue', lwd=2)

#plot scatterplot for experimental conditions
plot(log(filtered_tpm_table$TNZ_1_C/filtered_tpm_table$TNZ_1_T)
     , log(filtered_tpm_table$WS_2_C/filtered_tpm_table$WS_2_T)
     , xlab="log fold change of TNZ control and light", ylab = "log fold change of WS-2 control and light"
     #, main= "Comparison of log TPM between WS-2 control and light and TNZ control and light conditions"
     , pch=20, col=rgb(0.1,0.1,0.1,0.3))
abline(v=0,col='blue', lwd=2)
abline(h=0,col='blue', lwd=2)

###############################################################
######                  Complete TPM Table                #####
###############################################################

#read tpm table with ril lines and parents
tnz_ws2_ril_TPM_table <- read.csv("new_tpm_tables/tnz_ws2_ril_TPM_table_2.csv", row.names = 1)
#rownames(tnz_ws2_ril_TPM_table) <- tnz_ws2_ril_TPM_table[,1]
tnz_ws2_ril_TPM_table[,1] <- NULL

#function to calculate euclidean distance
euclidean <- function( b) {
  #print(b)
  sqrt(sum(((0-b[[1]])^2), ((0-b[[2]])^2)))}

#filter TPM table and remove rows where not all variables have at least a value of 5 and at least one variable in the row should be 20 or more
ril_filtered_tpm_table <- tnz_ws2_ril_TPM_table %>% filter_all(all_vars(.>=5)) %>% filter_all(any_vars(.>=20))

# calculate the log change for each of the ecotypes
filtered_tpm_table$condition_log_TPM_TNZ <- with(filtered_tpm_table, log(TNZ_1_C/TNZ_1_T))
filtered_tpm_table$condition_log_TPM_WS <- with(filtered_tpm_table,log(WS_2_C/WS_2_T))
#zip stuff to make calculating euclidean distance easier
zipped_list <- mapply(list, filtered_tpm_table$condition_log_TPM_TNZ, filtered_tpm_table$condition_log_TPM_WS, SIMPLIFY=F)

#calculate euclidean distance
filtered_tpm_table$condition_eucl_distance <- sapply(zipped_list, euclidean)



filtered_tpm_table$condition_log_TPM_TNZ_stat <- apply(filtered_tpm_table[5],1 ,function(x) if(x > 0) "+" else "-")
filtered_tpm_table$condition_log_TPM_WS_stat <- apply(filtered_tpm_table[6],1 ,function(x) if(x > 0) "+" else "-")

#data for experimental conditions
condition_quartile = c()


for (i in rownames(filtered_tpm_table)){
#divide the data in quadrants based on the scatterplot 
if (filtered_tpm_table[i, "condition_log_TPM_TNZ_stat"] == "-"){
  if (filtered_tpm_table[i,"condition_log_TPM_WS_stat"] == "-"){
    condition_quartile <- append(condition_quartile, "q3")
  }
  if (filtered_tpm_table[i,"condition_log_TPM_WS_stat"] == "+"){
    condition_quartile <- append(condition_quartile ,"q1")
  }
}

if (filtered_tpm_table [i,"condition_log_TPM_TNZ_stat"] == "+"){
  if (filtered_tpm_table[i,"condition_log_TPM_WS_stat"] == "-"){
    condition_quartile <- append(condition_quartile ,"q4")
  }
  if (filtered_tpm_table [i,"condition_log_TPM_WS_stat"] == "+"){
    condition_quartile <- append(condition_quartile ,"q2")
  }
}
  #end forloop
}

filtered_tpm_table$condition_quartile <- condition_quartile
new_filtered_tpm_table <- filtered_tpm_table[filtered_tpm_table$condition_eucl_distance>1,]
rownames_filtered_main_tpm <- rownames(new_filtered_tpm_table)


condition_df <- filtered_tpm_table[,5:10]
condition_df <- split(condition_df, condition_df$condition_quartile)
#write.csv(rownames_filtered_main_tpm,file="test_gprofiler.txt", quote = F, row.names=F)
condition_subdf <- ril_filtered_tpm_table[rownames(ril_filtered_tpm_table) %in% rownames_filtered_main_tpm,]
#sort data based on quartiles
condition_subdf_1 <- condition_subdf[rownames(condition_subdf) %in% rownames(condition_df[["q1"]]),]
condition_subdf_2 <- condition_subdf[rownames(condition_subdf) %in% rownames(condition_df[["q2"]]),]
condition_subdf_3 <- condition_subdf[rownames(condition_subdf) %in% rownames(condition_df[["q3"]]),]
condition_subdf_4 <- condition_subdf[rownames(condition_subdf) %in% rownames(condition_df[["q4"]]),]

#write quartile 2 and 3 to seperate txt files used for GO enrichment analysis
write.csv(rownames(condition_subdf_2),file="condition_subdf_2.txt", quote = F, row.names=F)
write.csv(rownames(condition_subdf_3),file="condition_subdf_3.txt", quote = F, row.names=F)

# rownames(condition_df[["q2"]])
# create heatmap with condtition data and look at if it is C or T
condition_names <- colnames(condition_subdf)
ct_condition <- c()
for (name in condition_names){
 ct_condition <- append( ct_condition, strsplit(name,'_')[[1]][[3]])
}
ct_condition_df<- data.frame(condition_names, ct_condition, row.names = 1)
pheatmap(condition_subdf, scale='row', show_colnames = F, show_rownames = F, 
         annotation_col = structure(ct_condition_df))

##################################################

#ril_filtered_tpm_table <- tnz_ws2_ril_TPM_table %>% filter_all(all_vars(.>=5)) %>% filter_all(any_vars(.>=20))
#calculate log tpm for light condition and control condition
filtered_tpm_table$samples_log_TPM_light <- with(filtered_tpm_table, log(WS_2_T/TNZ_1_T))
filtered_tpm_table$samples_log_TPM_control <- with(filtered_tpm_table,log(WS_2_C/TNZ_1_C))

# 
filtered_tpm_table$samples_log_TPM_light_stat <- apply(filtered_tpm_table[11],1 ,function(x) if(x > 0) "+" else "-")
filtered_tpm_table$samples_log_TPM_control_stat <- apply(filtered_tpm_table[12],1 ,function(x) if(x > 0) "+" else "-")

zipped_list <- mapply(list, filtered_tpm_table$samples_log_TPM_light, filtered_tpm_table$samples_log_TPM_control, SIMPLIFY=F)

#calculate euclidean distance
filtered_tpm_table$samples_eucl_distance <- sapply(zipped_list, euclidean)

samples_quartile <- c()
#filter data based on which quarter they were in in the scatterplot
for (i in rownames(filtered_tpm_table)){
  
  if (filtered_tpm_table[i, "samples_log_TPM_light_stat"] == "-"){
    if (filtered_tpm_table[i,"samples_log_TPM_control_stat"] == "-"){
      samples_quartile <- append(samples_quartile, "q3")
    }
    if (filtered_tpm_table[i,"samples_log_TPM_control_stat"] == "+"){
      samples_quartile <- append(samples_quartile ,"q4")
    }
  }
  
  if (filtered_tpm_table [i,"samples_log_TPM_light_stat"] == "+"){
    if (filtered_tpm_table[i,"samples_log_TPM_control_stat"] == "-"){
      samples_quartile <- append(samples_quartile ,"q1")
    }
    if (filtered_tpm_table [i,"samples_log_TPM_control_stat"] == "+"){
      samples_quartile <- append(samples_quartile ,"q2")
    }
  }
}

filtered_tpm_table$samples_quartile <- samples_quartile

new_filtered_tpm_table <- filtered_tpm_table[filtered_tpm_table$samples_eucl_distance>1,]

rownames_filtered_main_tpm <- rownames(new_filtered_tpm_table)


samples_df <- filtered_tpm_table[,11:16]
samples_df <- split(samples_df, samples_df$samples_quartile)
#write.csv(rownames_filtered_main_tpm,file="test_gprofiler.txt", quote = F, row.names=F)
samples_subdf <- ril_filtered_tpm_table[rownames(ril_filtered_tpm_table) %in% rownames_filtered_main_tpm,]
# split dataframe into quarter
samples_subdf_1 <- samples_subdf[rownames(samples_subdf) %in% rownames(samples_df[["q1"]]),]
samples_subdf_2 <- samples_subdf[rownames(samples_subdf) %in% rownames(samples_df[["q2"]]),]
samples_subdf_3 <- samples_subdf[rownames(samples_subdf) %in% rownames(samples_df[["q3"]]),]
samples_subdf_4 <- samples_subdf[rownames(samples_subdf) %in% rownames(samples_df[["q4"]]),]

selected_samples_for_elasticnet <- rownames(rbind(samples_subdf_1,samples_subdf_2,samples_subdf_3,samples_subdf_4))

selected_condition_for_elasticnet <- rownames(rbind(condition_subdf_1,condition_subdf_2,condition_subdf_3,condition_subdf_4))

combined_selected_samples <- c(selected_condition_for_elasticnet, selected_samples_for_elasticnet)

unique_combined_selected_samples <- unique(combined_selected_samples)

# these genes were used for elastic net
write.csv(unique_combined_selected_samples,file="genes_for_elasticnet.txt", quote = F, row.names=F)


# print quarter 2 and 3 gene names to txt file for gprofiler analysis
write.csv(rownames(samples_subdf_2),file="samples_subdf_2.txt", quote = F, row.names=F)
write.csv(rownames(samples_subdf_3),file="samples_subdf_3.txt", quote = F, row.names=F)

# create pheatmap
#subdf <- ril_filtered_tpm_table[rownames(ril_filtered_tpm_table) %in% rownames_filtered_main_tpm,]
sample_names <- colnames(samples_subdf)
ct_sample <- c()
for (name in sample_names){
  ct_sample <- append( ct_sample, strsplit(name,'_')[[1]][[3]])
}
ct_sample_df<- data.frame(sample_names, ct_sample, row.names = 1)
pheatmap(samples_subdf, scale='row', show_colnames = F, show_rownames = F, 
         annotation_col = structure(ct_sample_df))


