library(affy)

files = list.files("/Users/bencekover/Downloads/GSE31642_RAW/arrays", 
    full.names = TRUE)
affy.data = ReadAffy(filenames = files)
eset.mas5 = mas5(affy.data)


exprSet.nologs = exprs(eset.mas5)
# List the column (chip) names
colnames(exprSet.nologs)

colnames(exprSet.nologs) = c("wt1", "wt2", "fkh2d_1", "fkh2d_2", 
    "fkh2-S2A_1", "fkh2-S2A_2")


heatmap(exprSet.nologs)

#take logs
exprSet = log(exprSet.nologs, 2)
data.mas5calls = mas5calls(affy.data)
# Get the actual A/P calls
data.mas5calls.calls = exprs(data.mas5calls)


#calc means
wt_mean = apply(exprSet[, c("wt1","wt2")], 1, mean)
fkh2d_mean = apply(exprSet[, c("fkh2d_1", "fkh2d_2")], 1, mean)
fkh2s2a_mean = apply(exprSet[, c("fkh2-S2A_1", "fkh2-S2A_2")], 1, mean)

#calc log fold difs
wt_fkh2d = fkh2d_mean-wt_mean
wt_fkh2s2a = fkh2s2a_mean-wt_mean
fkh2d_fkh2s2a = fkh2s2a_mean-fkh2d_mean

#all data to include means and logfold changes
all.data = cbind(exprSet, wt_mean, fkh2d_mean, fkh2s2a_mean, wt_fkh2d, wt_fkh2s2a, fkh2d_fkh2s2a)
colnames(all.data)

#now calculating p values for all genes. welch test two sided
 wt_fkh2d_pval =  apply(exprSet, 1, function(x) {
    t.test(x[1:2], x[3:4], var.equal = FALSE)$p.value })
wt_fkh2s2a_pval =  apply(exprSet, 1, function(x) {
      t.test(x[1:2], x[5:6], var.equal = FALSE)$p.value })
fkh2d_fkh2s2a_pval =  apply(exprSet, 1, function(x) {
      t.test(x[3:4], x[5:6], var.equal = FALSE)$p.value })



# Concatenate all A/P calls for brain and liver
AP = apply(data.mas5calls.calls, 1, paste, collapse = "")

# Get the probsets where the 4 calls are not 'AAAA'
genes.present = names(AP[AP != "AAAAAA"])
# How many probetset/genes are present?
length(genes.present)

# Get all data for probesets that are present on at least on chip.
exprSet.present = exprSet[genes.present, ]

wt_fkh2d_pval_present =  wt_fkh2d_pval[genes.present]
wt_fkh2s2a_pval_present =  wt_fkh2s2a_pval[genes.present]
fkh2d_fkh2s2a_pval_present =  fkh2d_fkh2s2a_pval[genes.present]


wt_fkh2d_pval_present_adjusted = p.adjust(wt_fkh2d_pval_present, method = 'fdr')
wt_fkh2s2a_pval_present_adjusted = p.adjust(wt_fkh2s2a_pval_present, method = 'fdr')
fkh2d_fkh2s2a_pval_present_adjusted = p.adjust(fkh2d_fkh2s2a_pval_present, method = 'fdr')

#now order these

wt_fkh2d_pval_present_adjusted_ordered = wt_fkh2d_pval_present_adjusted[order(wt_fkh2d_pval_present_adjusted)]
wt_fkh2s2a_pval_present_adjusted_ordered = wt_fkh2s2a_pval_present_adjusted[order(wt_fkh2s2a_pval_present_adjusted)]
fkh2d_fkh2s2a_pval_present_adjusted_ordered = fkh2d_fkh2s2a_pval_present_adjusted[order(fkh2d_fkh2s2a_pval_present_adjusted)]


wt_fkh2s2a_pval_present_adjusted_ordered[1:10]








#rename the rownames based on this file /Users/bencekover/Downloads/TFS-Assets_LSG_Support-Files_Yeast_2-na36-annot-csv/Yeast_2.na36.annot.csv. heaer is 21th row
mapping_names = read.table("/Users/bencekover/Downloads/TFS-Assets_LSG_Support-Files_Yeast_2-na36-annot-csv/Yeast_2.na36.annot.csv", header = TRUE, sep = ",")
mapping_names = mapping_names[,c(1,7)]
colnames(mapping_names) = c("gene", "ID")
head(mapping_names)

#turn all.data to a data frame
all.data = as.data.frame(all.data)

#create new column
all.data_final= cbind(gene = rownames(all.data), all.data)

#iterate through gene column in all.data, find the row in mappingnames that matches in the first column, and replace the gene with  the value in column Transcript ID(Array Design)
for (i in 1:nrow(all.data_final)){
  all.data_final[i,1] = mapping_names[mapping_names$gene == all.data_final[i,1],2]
}
#save as csv
write.csv(all.data_final, file = "/Users/bencekover/Library/CloudStorage/OneDrive-UniversityCollegeLondon/MSci Bahler lab/S.-Pombe-biofilm/Bence folder/Analysis of microarray datasets/szilagyi_final.csv")

