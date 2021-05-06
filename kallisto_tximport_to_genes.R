library(rjson)
library(tximport)

#print out the system/R configuration
sessionInfo()
#get the arguments and create parameters from the json file
args <- commandArgs(trailingOnly = TRUE)
jObj <- fromJSON(file = args[1])
#summarize and print out the complete json object
summary(jObj)
jObj

#read in the design file.  the directory labels MUST have column header "sample"
samples <- read.table(jObj$design_file,header=TRUE)
samples
files <- file.path(jObj$baseDir,samples$sample,"kallisto_quant","abundance.h5")
files
names(files) <- samples$sample

tx2gene <- read.table(jObj$tx2gene_file)
summary(tx2gene)
txi <- tximport(files, type="kallisto", tx2gene=tx2gene)
summary(txi)

write.table(as.table(txi$counts),paste(jObj$output_file_prefix,'.gene.counts.txt',sep=""),sep='\t',col.names=NA,row.names=TRUE,quote=FALSE)
write.table(as.table(txi$abundance),paste(jObj$output_file_prefix,'.gene.tpm.txt',sep=""),sep='\t',col.names=NA,row.names=TRUE,quote=FALSE)

txit <- tximport(files, type="kallisto", txOut=TRUE)
summary(txit)
write.table(as.table(txit$counts),paste(jObj$output_file_prefix,'.transcript.counts.txt',sep=""),sep='\t',col.names=NA,row.names=TRUE,quote=FALSE)
write.table(as.table(txit$abundance),paste(jObj$output_file_prefix,'.transcript.tpm.txt',sep=""),sep='\t',col.names=NA,row.names=TRUE,quote=FALSE)
