rm(list = ls())
library(RMySQL)
con <- dbConnect(MySQL(), host="127.0.0.1", port=3306, user="root", password="11111111") 
dbSendQuery(con, "USE gse62944")
dbListTables(con)
setwd('G:\\GSE62944')
tumorCancerType2amples=read.table('GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt',sep = '\t',stringsAsFactors = F)
colnames(tumorCancerType2amples)=c('sampleID','CancerType')
tmp=lapply(unique(tumorCancerType2amples$CancerType), function(x){
  #x='PRAD';
  gene1="TP53";gene2="BRAC1";
  sqlTable=paste('tumor',x,'RPKM',sep='_')
  sqlQuery=paste0(' select * from ', sqlTable ,' where genesymbol = ',shQuote(gene1),' OR  genesymbol = ',shQuote(gene2))
  matrix2genes=dbGetQuery(con,sqlQuery)
  rownames(matrix2genes)=matrix2genes$geneSymbol
  matrix2genes=matrix2genes[,- match('geneSymbol',colnames(matrix2genes)) ]
  matrix2genes=t(matrix2genes)
  valueList1=as.numeric(matrix2genes[,gene1]);valueList2=as.numeric(matrix2genes[,gene2]);
  png(  paste0(gene1,'_and_',gene2,'_in_',x,'.SinalCor.png')  )
  plot(valueList1,valueList2,xlab=gene1,ylab=gene2)
  abline(lm(valueList2~valueList1),col='red')
  title(main =paste0("R2=",cor(valueList1,valueList2)))
  dev.off() 
  return(c(x,fivenum(valueList1),fivenum(valueList2),cor(valueList1,valueList2)))
})
write.csv(x = matrix(unlist(tmp),ncol=12,byrow = T),file = 'tumor.corration.csv')



normalCancerType2amples=read.table('GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt',sep = '\t',stringsAsFactors = F)
colnames(normalCancerType2amples)=c('sampleID','CancerType')

gene1="TP53";gene2="BRAC1"; 
sqlTable = 'normalrpkm';
sqlQuery=paste0(' select * from ', sqlTable ,' where genesymbol = ',shQuote(gene1),' OR  genesymbol = ',shQuote(gene2))
matrix2genes=dbGetQuery(con,sqlQuery)
rownames(matrix2genes)=matrix2genes$geneSymbol
matrix2genes=matrix2genes[,- match('geneSymbol',colnames(matrix2genes)) ]
matrix2genes=t(matrix2genes)
lapply(unique(normalCancerType2amples$CancerType), function(x){
  #x='THYM'; 
  sampleList=normalCancerType2amples[normalCancerType2amples$CancerType==x,1]
  sampleList=gsub("-",".", sampleList)
  matrix2genesType=matrix2genes[sampleList,]
  valueList1=as.numeric(matrix2genesType[,gene1]);valueList2=as.numeric(matrix2genesType[,gene2]);
  png(  paste0(gene1,'_and_',gene2,'_in_normal_',x,'.SinalCor.png')  )
  plot(valueList1,valueList2,xlab=gene1,ylab=gene2)
  abline(lm(valueList2~valueList1),col='red')
  title(main =paste0("R2=",cor(valueList1,valueList2)))
  dev.off() 
})


rm(list=ls())

searchGene = 'CBX6';
searchTable='tumor_brca_rpkm';

library(RMySQL)
con <- dbConnect(MySQL(), host="127.0.0.1", port=3306, user="root", password="11111111") 
dbSendQuery(con, "USE gse62944")
dbListTables(con)
query = paste0(' select * from ', searchTable ,' where genesymbol = ',shQuote(searchGene)) ;
expression_1=dbGetQuery(con,query)
expression_1=as.numeric(expression_1[,-1]);

query = paste0(' select geneSymbol from ', searchTable ) ;
allGenes=dbGetQuery(con,query)[,1]

cor_results <- matrix(unlist(lapply(allGenes, function(x){
  thisGene=x
  query = paste0(' select * from ', searchTable ,' where genesymbol = ',shQuote(thisGene)) ;
  expression_2=dbGetQuery(con,query)
  expression_2=as.numeric(expression_2[,-1]);
  tmp=cor.test(expression_1,expression_2);#str(tmp)
  return(c(thisGene,tmp$estimate,tmp$p.value))
}) ## end for lapply 
) ## end for unlist 
,ncol = 3, byrow =T) ## end for matrix




