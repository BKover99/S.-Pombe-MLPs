#install.packages("/Users/bencekover/Downloads/RFQTL.tar.gz",repos = NULL)

library(RFQTL)


library(randomForest)

library(qqman)
#loading in the genotype presence-absence matrix

#OLD
#d<-read.table(file = '/Users/bencekover/Library/CloudStorage/OneDrive-Personal/MSci Bahler lab/S.-Pombe-MLPs - Github/Bence folder/QTL mapping/SupplementaryDataset_S7_genotype.tsv', sep = '\t', header = TRUE)

#NEW
#d<-read.table(file = '/Users/bencekover/bioinformatics/projects/segregants/updated_genotype_matrix.tsv', sep = '\t', header = TRUE)
d<-read.table(file = '/Users/bencekover/bioinformatics/projects/segregants/updated_genotype_matrix_final_pos.tsv', sep = '\t', header = TRUE)

d<-d[c(1:(nrow(d)-12)),]
genotype1<-d[,-c(1:4)]
genotype1<-data.matrix(genotype1)
genotype1<-t(genotype1)
#loading in the filtering measurements
phenotype1<-read.csv("/Users/bencekover/Library/CloudStorage/OneDrive-Personal/MSci Bahler lab/S.-Pombe-MLPs - Github/Bence folder/Image processing - Adhesion assay/segregants_emm_p.csv")

phenotype1<-phenotype1[which(phenotype1$strain %in% rownames(genotype1)),]
strainNames1 <- phenotype1$strain
phenotype1 <- phenotype1$ratio
sampleInfo1 <- sapply(strainNames1,FUN=function(x){
  which(rownames(genotype1)==x)
})
mode(genotype1) <- "integer"
mappingData1 <- preMap(genotype=genotype1,
                       phenotype=phenotype1,
                       sampleInfo=sampleInfo1,
                       scale=T,
                       )


#Get accurate realscores by calculating them 20x and averaging them

r=rfMapper(mappingData = mappingData1,
           permute = F,
           nforest = 100,#was 100
           ntree = 100)

#Get accurate p-values by doing 10000 permutations (100x35)
#Here the wd is a folder in which your permutations will go

#To test old
setwd("/Users/bencekover/Library/CloudStorage/OneDrive-Personal/MSci Bahler lab/S.-Pombe-MLPs - Github/Bence folder/QTL mapping/permutations_emm_p/permutations")

permutedScores1 <- rfMapper(mappingData = mappingData1,
                            permute = T,
                            nforest = 100,#was 100
                            ntree = 100,
                            nPermutations=2000, #was 100
                            file="filter_permutt.RData",
                            nCl=4,
                            clType="SOCK")

#Below the wd is the folder in which the permutations folder is (not IN the permutation folder)
setwd("/Users/bencekover/Library/CloudStorage/OneDrive-Personal/MSci Bahler lab/S.-Pombe-MLPs - Github/Bence folder/QTL mapping/permutations_emm_p/")
#The path below is the permutations folder

pValues1 <- pEst(path="permutations/",
                 scores=r,
                 markersPerIteration = 350,
                 printProg = T,
                 pCorrection = "none")
pValuesX1 <- pValues1[mappingData1$genotype2group]




#df for export to make Manhattan plot in Python
df <- data.frame (haplogroup  = mappingData1$genotype2group,
                  chr = d$chromosome,
                  p_val= pValues1[mappingData1$genotype2group])

write.csv(df,"/Users/bencekover/Library/CloudStorage/OneDrive-Personal/MSci Bahler lab/S.-Pombe-MLPs - Github/Bence folder/Manhattan plot/dataresults_emm_p.csv")

#These steps (until markerPositions) are not useful for the way I presented the results
chrVec1 <- d$chromosome


QTL_list1 <- QTLgrouper(pmat = pValuesX1,
                        sigThreshold = 0.001,
                        corThreshold = 0.99,
                        distThreshold = 1,
                        genotype = genotype1,
                        chrVec = chrVec1)

#Plotting the marker positions

markerPositions1 <- d[,c(1,2)]
markerPositions1$chromosome<-gsub("chromosome_1",1,markerPositions1$chromosome)
markerPositions1$chromosome<-gsub("chromosome_2",2,markerPositions1$chromosome)
markerPositions1$chromosome<-gsub("chromosome_3",3,markerPositions1$chromosome)

markerPositions1[,3]=markerPositions1$position


#the writeQTL is what is in the tutorial but I didn't find it useful as a visualisation of results
writeQTL(QTLlist = QTL_list1,traitNames = "Flocc",markerPositions = markerPositions1,path="myResults_new_filtering.qtl")
qtl1 <- readQTL(path = "myResults_new_filtering.qtl")
qtl1





#Presenting results
#barplot
barplot(-log10(pValuesX1))
abline(h=-log10(0.05/length(pValues1)),col="red")


#manhattan plot
results.tab=markerPositions1[,-3]
results.tab[,3]=pValuesX1
colnames(results.tab)=c("CHR","BP","P")

results.tab=results.tab[-which(is.na(results.tab)),]
results.tab$CHR=as.numeric(unlist(results.tab$CHR))
results.tab[,4]=rep("snp",length(results.tab$CHR))
colnames(results.tab)[4]="SNP"
manhattan(results.tab,ylim = c(0,7), suggestiveline = F, genomewideline = F,cex = 0.7)
abline(h=-log10(0.05/length(pValues1)),col="red")
abline(h=-log10(1/20000),col="blue")
text(x=200,y=-log10(0.05/length(pValues1))+0.15,labels="                   p=7.2e-05",col="red",cex=0.7)

#from the following, determine the regions to use for reg_map below
snps = markerPositions1[which(-log10(pValuesX1)>=3),] 
results = d[row.names(snps),]
write.csv(results,"results_emm_p.csv")
plot(genotype1[,qtl1[[1]]$mostSignificantPredictor], phenotype1)
plot(genotype1[,qtl1[[4]]$mostSignificantPredictor], phenotype1)


d[qtl1[[4]]$predictors,]

#manually
#mapping regions on chromosomes
chromosomes=c("chr1","chr2","chr3")
start=c(1,1,1)
end=c(5579133,4539804,2452883)
chr_map=matrix(c(chromosomes, start, end), ncol=3, byrow=F)
chr_map=data.frame(chr_map)
regions=c("1")
chr_names=c("chr2")
reg_start=c(2528629)
reg_end=c(2536217)
reg_map=matrix(c(regions,chr_names, reg_start, reg_end), ncol=4, byrow=F)
reg_map=data.frame(reg_map)
write.table(chr_map,"chr_map.txt",sep="\t",col.names=F,row.names=F)
write.table(reg_map,"reg_map.txt",sep="\t",col.names=F,row.names=F)
library(chromoMap)
chromoMap("chr_map.txt","reg_map.txt",segment_annotation=T)
