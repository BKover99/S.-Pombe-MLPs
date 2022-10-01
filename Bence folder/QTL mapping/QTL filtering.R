d<-read.table(file = 'SupplementaryDataset_S7_genotype.tsv', sep = '\t', header = TRUE)
library(RFQTL)
genotype2<-d[,-c(1:4)]
genotype2<-data.matrix(genotype2)
genotype2<-t(genotype2)
phenotype2<-read.csv('Biofilm_image_processing/Biofilm_image_processing/filtering_phenotypes.csv')
phenotype2<-phenotype2[,-1]
phenotype2<-phenotype2[which(phenotype2$Strain %in% rownames(genotype2)),]
strainNames2 <- phenotype2$Strain
phenotype2 <- phenotype2$X.flocc
sampleInfo2 <- sapply(strainNames2,FUN=function(x){
  which(rownames(genotype2)==x)
})
mappingData2 <- preMap(genotype=genotype2,
                       phenotype=phenotype2,
                       sampleInfo=sampleInfo2,
                       scale=T)
library(randomForest)
realScores2 <- rfMapper(mappingData = mappingData2,
                        permute = F,
                        nforest = 100,
                        ntree = 150)
permutedScores2 <- rfMapper(mappingData = mappingData2,
                            permute = T,
                            nforest = 100,
                            ntree = 150,
                            nPermutations=10,
                            file="filt_permut1.RData",
                            nCl=4,
                            clType="SOCK")
pValues2 <- pEst(path="filtering_permutations/",
                 scores=realScores2,
                 markersPerIteration = 350,
                 printProg = T,
                 pCorrection = "none")
pValuesX2 <- pValues2[mappingData2$genotype2group]
chrVec2 <- d$chromosome
QTL_list2 <- QTLgrouper(pmat = pValuesX2,
                        sigThreshold = 0.1,
                        corThreshold = 0.8,
                        distThreshold = 9,
                        genotype = genotype2,
                        chrVec = chrVec2)
markerPositions2 <- d[,c(1,2)]
markerPositions2$chromosome<-gsub("chromosome_1",1,markerPositions2$chromosome)
markerPositions2$chromosome<-gsub("chromosome_2",2,markerPositions2$chromosome)
markerPositions2$chromosome<-gsub("chromosome_3",3,markerPositions2$chromosome)
markerPositions2[,3]=markerPositions2$position
writeQTL(QTLlist = QTL_list2,traitNames = "Flocc",markerPositions = markerPositions2,path="myResults2.qtl")
qtl2 <- readQTL(path = "myResults2.qtl")
qtl2