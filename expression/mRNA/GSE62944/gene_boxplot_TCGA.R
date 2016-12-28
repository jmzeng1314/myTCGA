rm(list=ls())

searchGene = 'VCX3B';
searchTable1='tumor_gbm_rpkm';
searchTable2='tumor_lgg_rpkm';

library(RMySQL)
con <- dbConnect(MySQL(), host="127.0.0.1", port=3306, user="root", password="11111111") 
dbSendQuery(con, "USE gse62944")
dbListTables(con)
query = paste0(' select * from ', searchTable1 ,' where genesymbol = ',shQuote(searchGene)) ;
gbm=dbGetQuery(con,query)
query = paste0(' select * from ', searchTable2 ,' where genesymbol = ',shQuote(searchGene)) ;
lgg=dbGetQuery(con,query)

gbm=as.numeric(gbm[,-1]);gbm=data.frame(value=gbm,type='gbm')
lgg=as.numeric(lgg[,-1]);lgg=data.frame(value=lgg,type='lgg')
dat1= rbind(gbm,lgg)

boxplot( value ~  type, data = dat1, lwd = 2, ylab = 'value')
stripchart(value ~ type, vertical = TRUE, data = dat1, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

 
 
sqlTable = 'normalrpkm';
sqlQuery=paste0(' select * from ', sqlTable ,' where genesymbol = ',shQuote(searchGene))
normalExpression=dbGetQuery(con,sqlQuery)
normalExpression= normalExpression[,-length(normalExpression)] 
normalExpression = data.frame(sampleID=names(normalExpression),
                              values=as.numeric(normalExpression)
                              )
normalCancerType2amples=dbGetQuery(con,'select * from normalcancertype2amples')
normalCancerType2amples$sampleID=gsub("-",".", normalCancerType2amples$sampleID)
dat2 = merge(normalExpression,normalCancerType2amples,by='sampleID')

boxplot( values ~  CancerType, data = dat2, lwd = 2, ylab = 'values',las=2,main=searchGene)
stripchart(values ~ CancerType, vertical = TRUE, data = dat2, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
