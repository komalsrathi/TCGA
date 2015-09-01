################################ prep up the TCGA annotation ##############################################                                                 
system('grep \'gene\' TCGA.Sept2010.09202010.gaf | cut -f16,17 > TCGA.Sept2010.09202010.gene.gaf') # [Gene & GeneLocus]
tcga = read.delim('TCGA.Sept2010.09202010.gene.gaf',header=F)
s <- strsplit(as.character(tcga$V2), ';')
final = data.frame(V1 = rep(tcga$V1, sapply(s, length)), V2 = unlist(s))
final = cbind(final,colsplit(final$V2,pattern=":",names=c("chr","coord","strand")))
final$V1 = sub("[|][0-9]of[0-9]","",final$V1)
final = unique(final)
final = cbind(final,colsplit(final$coord,pattern="-",names=c("start","end")))
final = transform(final, start = ifelse(strand == '+', start, end), end = ifelse(strand == '+', end, start))
final = final[,c(1,3,5,6,7)]
final = ddply(final, .(V1,chr,strand), summarise, start= min(start), end= max(end))
final$score = 0
final = final[,c(2,4,5,1,6,3)]
write.table(final,'TCGA.Sept2010.09202010.gene.bed',row.names=F,col.names=T,sep="\t",quote=F)
################################ prep up the TCGA annotation ##############################################   
