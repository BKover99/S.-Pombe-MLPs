setwd("/Users/celestecohen/Downloads")
d<-read.table(file = 'SupplementaryDataset_S7_genotype.tsv', sep = '\t', header = TRUE)
library(RFQTL)
genotype1<-d[,-c(1:4)]
genotype1<-data.matrix(genotype1)
genotype1<-t(genotype1)
phenotype1<-read.csv('Biofilm_bioinformatics/Biofilm_image_processing/washing_phenotypes.csv')
phenotype1<-phenotype1[,-1]
phenotype1<-phenotype1[which(phenotype1$strain %in% rownames(genotype1)),]
strainNames1 <- phenotype1$strain
phenotype1 <- phenotype1$ratio
sampleInfo1 <- sapply(strainNames1,FUN=function(x){
  which(rownames(genotype1)==x)
})
mappingData1 <- preMap(genotype=genotype1,
                      phenotype=phenotype1,
                      sampleInfo=sampleInfo1,
                      scale=T)
library(randomForest)

#Get accurate realscores by calculating them 20x and averaging them
RS1=c()
for(i in 1:20){
  r=rfMapper(mappingData = mappingData1,
             permute = F,
             nforest = 100,
             ntree = 150)
  RS1=c(RS1,r)
}
RS1_df=matrix(RS1, ncol=20, byrow=F)
RS1_df=data.frame(RS1_df)
realScores1<-rowMeans(RS1_df)

#Get accurate p-values by doing 3500 permutations (100x35)
#Here the wd is a folder in which your permutations will go
setwd("/Users/celestecohen/Downloads/Biofilm_bioinformatics/washing_permutations/permutations")
for(i in 1:35){
  permutedScores1 <- rfMapper(mappingData = mappingData1,
                              permute = T,
                              nforest = 100,
                              ntree = 150,
                              nPermutations=100,
                              file=paste("wash_permut",i,".RData",sep=""),
                              nCl=4,
                              clType="SOCK")
}

#Below the wd is the folder in which the permutations folder is (not IN the permutation folder)
setwd("/Users/celestecohen/Downloads/Biofilm_bioinformatics")
#The path below is the permutations folder
pValues1 <- pEst(path="washing_permutations/permutations/",
                scores=realScores1,
                markersPerIteration = 350,
                printProg = T,
                pCorrection = "none")
pValuesX1 <- pValues1[mappingData1$genotype2group]

#These steps (until markerPositions) are not useful for the way I presented the results
chrVec1 <- d$chromosome
QTL_list1 <- QTLgrouper(pmat = pValuesX1,
                       sigThreshold = 0.1,
                       corThreshold = 0.8,
                       distThreshold = 9,
                       genotype = genotype1,
                       chrVec = chrVec1)

#this is useful
markerPositions1 <- d[,c(1,2)]
markerPositions1$chromosome<-gsub("chromosome_1",1,markerPositions1$chromosome)
markerPositions1$chromosome<-gsub("chromosome_2",2,markerPositions1$chromosome)
markerPositions1$chromosome<-gsub("chromosome_3",3,markerPositions1$chromosome)
markerPositions1[,3]=markerPositions1$position

#the writeQTL is what is in the tutorial but I didn't find it useful as a visualisation of results
writeQTL(QTLlist = QTL_list1,traitNames = "Flocc",markerPositions = markerPositions1,path="myResults1.qtl")
qtl1 <- readQTL(path = "myResults1.qtl")
qtl1


#Presenting results
#barplot
barplot(-log10(pValuesX1))
abline(h=-log10(0.01),col="red")

#manhattan plot
results.tab=markerPositions1[,-3]
results.tab[,3]=pValuesX1
colnames(results.tab)=c("CHR","BP","P")
library(qqman)
results.tab=results.tab[-which(is.na(results.tab)),]
results.tab$CHR=as.numeric(unlist(results.tab$CHR))
results.tab[,4]=rep("snp",length(results.tab$CHR))
colnames(results.tab)[4]="SNP"
manhattan(results.tab,ylim = c(0,4), suggestiveline = F, genomewideline = F,cex = 0.7)
abline(h=-log10(0.01),col="red")
text(x=1,y=-log10(0.01)+0.1,labels="       p=0.01",col="red",cex=0.7)

#from the following, determine the regions to use for reg_map below
markerPositions1[which(-log10(pValuesX1)>=2),] 


#mapping regions on chromosomes
chromosomes=c("chr1","chr2","chr3")
start=c(1,1,1)
end=c(5579133,4539804,2452883)
chr_map=matrix(c(chromosomes, start, end), ncol=3, byrow=F)
chr_map=data.frame(chr_map)
regions=c("1","2","3","4","5")
chr_names=c("chr1","chr2","chr2","chr2","chr3")
reg_start=c(2055773,2202682,2334988,2516209,210341)
reg_end=c(2055773,2202682,2363654,2640722,210364)
reg_map=matrix(c(regions,chr_names, reg_start, reg_end), ncol=4, byrow=F)
reg_map=data.frame(reg_map)
write.table(chr_map,"chr_map.txt",sep="\t",col.names=F,row.names=F)
write.table(reg_map,"reg_map.txt",sep="\t",col.names=F,row.names=F)
library(chromoMap)
chromoMap("chr_map.txt","reg_map.txt",segment_annotation=T)

