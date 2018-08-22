#!/usr/local/bin/Rscript --verbose

setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants")

###setting things up ----
#library(vegan)
#library(gtools)
library(dada2)
#library(phyloseq)

#the script will run for either the fungi OR bacteria
#(note that I modified the name of the fungi file to follow the same nomenclature as the bacteria files (done using "rename  's/9_P/9_Fun_P/' *fastq.gz", "rename  's/9_R/9_Fun_R/' *fastq.gz" AND "rename  's/9_T/9_Fun_T/' *fastq.gz")

#list parameters 
reference.taxonomy = "reference_material/silva_nr_v132_train_set.fa"
path = "sequences/bacteria/"
group = "bacteria"
#r1_12= c(195,203,211,219,227,234,242,250,258,266,274,281)
time=Sys.time();names(time) = "start"

#####list fastq files -----
system(paste("ls -1 ",path,"raw/*_R1.fastq.gz >results/R1.list.fastq.gz",sep = ""))
list.fastq = read.table("results/R1.list.fastq.gz",stringsAsFactors = F);colnames(list.fastq) = "R1zip"

#order by R_,P_,T_
R_ = list.fastq[regexpr("Bac_R",list.fastq[,1])>0,1]; P_ = list.fastq[regexpr("Bac_P",list.fastq[,1])>0,1]; T_ = list.fastq[regexpr("Bac_T",list.fastq[,1])>0,1]
list.fastq[,1] = c(R_,P_,T_)

#gsub the names
list.fastq$R2zip = gsub("_R1.","_R2.",list.fastq$R1zip,fixed = TRUE)
list.fastq$R1filt = gsub(".gz",".filt.gz",list.fastq$R1zip,fixed = TRUE);list.fastq$R1filt = gsub("raw","filtered",list.fastq$R1filt)
list.fastq$R2filt = gsub(".gz",".filt.gz",list.fastq$R2zip,fixed = TRUE);list.fastq$R2filt = gsub("raw","filtered",list.fastq$R2filt)
sample.names = unlist(strsplit(list.fastq[,1],split = "_"))[seq(4,nrow(list.fastq)*5,by = 5)]

r1_96 = c(1:length(sample.names))[as.numeric(gsub("R","",sample.names))<49][1:48] #tomato roots.
nb_seq = 1:288 #nb of sequences to process. #

print(nb_seq)

time=c(time,Sys.time()); names(time)[length(time)] = "list"

#mkdir if necessary
if(file.exists(paste(path,"filtered/",sep = "")) == F) system(paste("mkdir ",path,"filtered/",sep = ""))


#filter and trim
out <- filterAndTrim(list.fastq$R1zip[nb_seq], list.fastq$R1filt[nb_seq], list.fastq$R2zip[nb_seq], list.fastq$R2filt[nb_seq],
                     truncQ=3,trimLeft=c(0,0),maxEE = c(10,10),minLen = 100,multithread=T,verbose =T)

print(paste("Done filter&trim",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "filter&trim"

#plot qual before and after
#plotQualityProfile(list.fastq$R1zip[1:4],aggregate =T)
#plotQualityProfile(list.fastq$R1filt[1:4],aggregate =T)

#error rates of sequences
error_R1 <- learnErrors(list.fastq$R1filt[c(1:4,97:100,193:196)],multithread = TRUE)
error_R2 <- learnErrors(list.fastq$R2filt[c(1:4,97:100,193:196)],multithread = TRUE)
#error_R1 <- learnErrors(list.fastq$R1filt[1:4],multithread = TRUE)
#error_R2 <- learnErrors(list.fastq$R2filt[1:4],multithread = TRUE)
#plotErrors(error_R1, nominalQ=TRUE)

print(paste("Done error",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "error"

#dereplication
derep_R1 <- derepFastq(list.fastq$R1filt[nb_seq])
derep_R2 <- derepFastq(list.fastq$R2filt[nb_seq])
names(derep_R1) <- sample.names[nb_seq]
names(derep_R2) <- sample.names[nb_seq]

print(paste("Done derep",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "derep"

#Sample Inference (reverse, R1 and forward, R2)
dada_R1 <- dada(derep_R1,err = error_R1,pool=FALSE,multithread = TRUE)
dada_R2 <- dada(derep_R2,err = error_R2,pool=FALSE,multithread = TRUE)

print(paste("Done dada inference",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "sample_inference"

#merging of PE
merged <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2,minOverlap = 20, maxMismatch = 3)

print(paste("Done merging",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "merging"

#Identify ASVs
asv <- makeSequenceTable(merged)

print(paste("Done ID ASVs",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "asv"

#Remove chimeras
asv.nochim <- removeBimeraDenovo(asv, method = "pooled", multithread = TRUE,verbose = TRUE) 
asv.nochim.bin <- ifelse(asv.nochim>0,1,0) 

asv.nochim.simplified = asv.nochim
colnames(asv.nochim.simplified) = paste("ASV",c(1:ncol(asv.nochim)),sep = "")

print(paste("Done remove chimeras",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "chimeras"

#summary stats 
getN <- function(x) sum(getUniques(x))
summary.stats <- cbind(out, # input & filtered
               round(((out[,2]/out[,1])*100),2), # % (filtered / input) 
               sapply(merged, getN), # Merged 
               round(((sapply(merged, getN)/out[,1])*100),2), # % (Merged / Input)
               rowSums(asv.nochim),# Non-chimeric
               round(((rowSums(asv.nochim)/out[,1])*100),2),# % (Non-chimeric / Input)
               rowSums(asv.nochim.bin)) # Number of ASVs per sample 

colnames(summary.stats) <- c("Input","Filter","%Filt/In","Merge","%Mer/In","Non-chim","%Non-chim/In","#ASV/sample") # Column names
rownames(summary.stats) <- sample.names[nb_seq] # Row names
summary.stats<-as.data.frame(summary.stats)
print(summary.stats)

time=c(time,Sys.time()); names(time)[length(time)] = "summary_stats"

#assign taxonomy
asv.taxo <- assignTaxonomy(asv.nochim,refFasta = reference.taxonomy,minBoot = 50, multithread=TRUE)
asv.taxo = cbind(asv.taxo,rownames(asv.taxo))
rownames(asv.taxo) = paste("ASV",c(1:nrow(asv.taxo)),sep = "")

print(paste("Done Assign taxo",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "taxo"

#times 
timed = data.frame(round(time[2:length(time)] - time[1:(length(time)-1)],2))
colnames(timed) = "Time"
timed = rbind(timed,sum(timed[,1]));row.names(timed)[11] = "Total time"
print(timed)

#save results ----
write.table(asv.nochim.simplified,paste("results/asv.nochim.",group,sep = ""))
write.table(asv.taxo,paste("results/asv.taxo.",group,sep = ""))


