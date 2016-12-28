rm(list = ls())
library(RMySQL)
con <- dbConnect(MySQL(), host="127.0.0.1", port=3306, user="root", password="11111111") 
dbSendQuery(con, "USE gse62944")
dbListTables(con)
setwd('G:\\GSE62944')
tumorCancerType2amples=read.table('GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt',sep = '\t',stringsAsFactors = F)
colnames(tumorCancerType2amples)=c('sampleID','CancerType')
dbWriteTable(con, 'tumorCancerType2amples', tumorCancerType2amples, append=F,row.names=F) 
normalCancerType2amples=read.table('GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt',sep = '\t',stringsAsFactors = F)
colnames(normalCancerType2amples)=c('sampleID','CancerType')
dbWriteTable(con, 'normalCancerType2amples', normalCancerType2amples, append=F,row.names=F)

Clinical_Variables=read.table('GSE62944_06_01_15_TCGA_24_548_Clinical_Variables_9264_Samples.txt',sep = '\t',stringsAsFactors = F)

normalRPKM=read.table('GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FPKM.txt.gz',sep = '\t',stringsAsFactors = F,header = T)
colnames(normalRPKM)[1]='geneSymbol'
rownames(normalRPKM)=normalRPKM$geneSymbol
normalRPKM=normalRPKM[,-1]
normalRPKM=round( as.matrix(normalRPKM),3)  
normalRPKM=as.data.frame(normalRPKM)
normalRPKM$geneSymbol = rownames(normalRPKM)
dbWriteTable(con, 'normalRPKM', normalRPKM, append=F,row.names=F)


tumorRPKM=read.table('GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt.gz',sep = '\t',stringsAsFactors = F,header = T)
colnames(tumorRPKM)[1]='geneSymbol'
rownames(tumorRPKM)=tumorRPKM$geneSymbol
tumorRPKM=tumorRPKM[,-1]
tumorRPKM=round( as.matrix(tumorRPKM),3)  
tumorRPKM=as.data.frame(tumorRPKM)
tumorRPKM$geneSymbol = rownames(tumorRPKM)
#load(file = 'tumorRPKM.rData')
lapply(unique((tumorCancerType2amples$CancerType)), function(x){
  #x='PRAD';
  sampleList=tumorCancerType2amples[tumorCancerType2amples$CancerType==x,1]
  sampleList=gsub("-",".", sampleList)
  tmpMatrix=tumorRPKM[,c('geneSymbol',sampleList)]
  dbWriteTable(con, paste('tumor',x,'RPKM',sep='_'), tmpMatrix, append=F,row.names=F)
  
})
