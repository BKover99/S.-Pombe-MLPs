library(affy)

files = list.files("/Users/bencekover/Library/CloudStorage/OneDrive-UniversityCollegeLondon/MSci Bahler lab/S.-Pombe-biofilm/external data/Li et al 2013/GSE43248_RAW", 
    full.names = TRUE)
affy.data = ReadAffy(filenames = files)
eset.mas5 = mas5(affy.data)


exprSet.nologs = exprs(eset.mas5)
# List the column (chip) names
colnames(exprSet.nologs)

colnames(exprSet.nologs) = c("WT1", "WT2", "WT3", "rpl3201_1", "rpl3201_2", "rpl3201_3", 
    "rpl3202_1", "rpl3202_2", "rpl3202_3")
 

heatmap(exprSet.nologs)
  
#take logs
exprSet = log(exprSet.nologs, 2)
data.mas5calls = mas5calls(affy.data)
# Get the actual A/P calls
data.mas5calls.calls = exprs(data.mas5calls)


#calc means
wt_mean = apply(exprSet[, c("WT1", "WT2", "WT3")], 1, mean)
rpl3201_mean = apply(exprSet[, c("rpl3201_1", "rpl3201_2", "rpl3201_3")], 1, mean)
rpl3202_mean = apply(exprSet[, c("rpl3202_1", "rpl3202_2", "rpl3202_3")], 1, mean)

#calc log fold difs
wt_rpl3201 = rpl3201_mean - wt_mean
wt_rpl3202 = rpl3202_mean - wt_mean 
rpl3201_rpl3202 =  rpl3202_mean - rpl3201_mean

#all data to include means and logfold changes
all.data = cbind(exprSet, wt_mean, rpl3201_mean, rpl3202_mean, wt_rpl3201, wt_rpl3202, rpl3201_rpl3202)
colnames(all.data)

#now calculating p values for all genes. welch test two sided
#iterate through table and do t test between columns
wt_rpl3201_pval = c()
for (i in 1:length(exprSet[,1])) {
  wt_rpl3201_pval = c(wt_rpl3201_pval, t.test(exprSet[i,1:3], exprSet[i,4:6], var.equal = FALSE)$p.value)
  #add name
    
}
names(wt_rpl3201_pval) = rownames(exprSet)
#show smallest p values
sort(wt_rpl3201_pval)[1:10]
#which gene
names(sort(wt_rpl3201_pval)[1:10])

# Concatenate all A/P calls for brain and liver
AP = apply(data.mas5calls.calls, 1, paste, collapse = "")

# Get the probsets where the 4 calls are not 'AAAA'
genes.present = names(AP[AP != "AAAAAA"])
# How many probetset/genes are present?
length(genes.present)

# Get all data for probesets that are present on at least on chip.
exprSet.present = exprSet[genes.present, ]

wt_rpl3201_pval_present =  wt_rpl3201_pval[genes.present]
wt_rpl3202_pval_present =  wt_rpl3202_pval[genes.present]
rpl3201_rpl3202_pval_present =  rpl3201_rpl3202[genes.present]


wt_rpl3201_pval_present_adjusted= p.adjust(wt_rpl3201_pval_present, method = 'fdr')

#show hist of p values on log scale on 100 bins

hist(log10(wt_rpl3201_pval_present), breaks = 100, col = "red", main = "p values for wt vs rpl3201", xlab = "p value")
wt_rpl3202_pval_present_adjusted= p.adjust(wt_rpl3202_pval_present, method = 'fdr')
rpl3201_rpl3202_pval_present_adjusted= p.adjust(rpl3201_rpl3202_pval_present, method = 'fdr')

#now order these

wt_rpl3201_pval_present_adjusted_ordered = wt_rpl3201_pval_present[order(wt_rpl3201_pval_present)]
wt_rpl3202_pval_present_adjusted_ordered = wt_rpl3202_pval_present[order(wt_rpl3202_pval_present)]
rpl3201_rpl3202_pval_present_adjusted_ordered = rpl3201_rpl3202_pval_present[order(rpl3201_rpl3202_pval_present)]


wt_rpl3201_pval_present



hist(wt_rpl3201_pval_present_adjusted, breaks = 100, col = "red", main = "p values for wt vs rpl3201", xlab = "p value")


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
write.csv(all.data_final, file = "/Users/bencekover/Library/CloudStorage/OneDrive-UniversityCollegeLondon/MSci Bahler lab/S.-Pombe-biofilm/Bence folder/Analysis of microarray datasets/li_final.csv")


all.data
mapping_names
