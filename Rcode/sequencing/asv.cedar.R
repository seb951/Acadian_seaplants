#!//cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/r/3.5.0/bin/Rscript --verbose

###setting things up ----
args = commandArgs(TRUE)
group = args[1]
nb_seq = seq(as.numeric(strsplit(args[2],":")[[1]][1]),as.numeric(strsplit(args[2],":")[[1]][2]),by = 1)
if(is.na(group)) {print("you didn't specify a proper group (e.g. bacteria_root) to analyse");quit(save="no")} else print(paste("analysing group:",group))
if(is.na(nb_seq[1])) {print("you didn't specify a proper number of sequences (eg.1:4) to analyse");quit(save="no")} else {print(paste("analysing sequences:"));print(nb_seq)}
  
setwd("/home/renaut/scratch/Acadian_seaplants")

library(dada2)
#library(phyloseq)

#the script will run for either the fungi OR bacteria
#(note that I modified the name of the fungi file to follow the same nomenclature as the bacteria files (done using "rename  's/9_P/9_Fun_P/' *fastq.gz", "rename  's/9_R/9_Fun_R/' *fastq.gz" AND "rename  's/9_T/9_Fun_T/' *fastq.gz")
#(note that I ran these two commands on the raw fastq files as the names were oddly formulated in some samples because of the command above that fucked things up....):
#:eg.: MI.M03992_0111.001.FLD0103.SCI016829_Fun_T79_Fun_R1.fastq.gz
#rename s/_Fun_R1.fastq.gz/_R1.fastq.gz/ *_Fun_R1.fastq.gz OR  #rename _Fun_R1.fastq.gz _R1.fastq.gz  *_Fun_R1.fastq.gz on Linux computer cluster.
#rename s/_Fun_R1.fastq.gz/_R1.fastq.gz/ *_Fun_R1.fastq.gz

#list parameters 
reference.taxonomy = ifelse(strsplit(group,"_")[[1]][1] == "bacteria", "reference_material/silva_nr_v132_train_set.fa","reference_material/unite_singleton/sh_general_release_dynamic_s_01.12.2017_reformated.fasta")
path = ifelse(strsplit(group,"_")[[1]][1] == "bacteria","sequences/bacteria/","sequences/fungi/")

print(paste("reference taxo is:",reference.taxonomy))
print(paste("path to sequence is:",path))

time=Sys.time();names(time) = "start"
print(paste("Starting time is:",time))


#####list fastq files -----
system(paste("ls -1 ",path,"raw/*_R1.fastq.gz >results/R1.list.fastq.gz",sep = ""))
list.fastq = read.table("results/R1.list.fastq.gz",stringsAsFactors = F);colnames(list.fastq) = "R1zip"

#order by R_,P_,T_
group_names = ifelse(strsplit(group,"_")[[1]][1] == "bacteria","Bac","Fun")
R_ = list.fastq[regexpr(paste(group_names,"_R",sep = ""),list.fastq[,1])>0,1]; P_ = list.fastq[regexpr(paste(group_names,"_P",sep = ""),list.fastq[,1])>0,1]; T_ = list.fastq[regexpr(paste(group_names,"_T",sep = ""),list.fastq[,1])>0,1]
list.fastq[,1] = c(R_,P_,T_)

#gsub the names
list.fastq$R2zip = gsub("_R1.","_R2.",list.fastq$R1zip,fixed = TRUE)
list.fastq$R1filt = gsub(".gz",".filt.gz",list.fastq$R1zip,fixed = TRUE);list.fastq$R1filt = gsub("raw","filtered",list.fastq$R1filt)
list.fastq$R2filt = gsub(".gz",".filt.gz",list.fastq$R2zip,fixed = TRUE);list.fastq$R2filt = gsub("raw","filtered",list.fastq$R2filt)
sample.names = unlist(strsplit(list.fastq[,1],split = "_"))[seq(4,nrow(list.fastq)*5,by = 5)]

time=c(time,Sys.time()); names(time)[length(time)] = "list"

#mkdir if necessary
if(file.exists(paste(path,"filtered/",sep = "")) == F) system(paste("mkdir ",path,"filtered/",sep = ""))

#filter and trim
out <- filterAndTrim(list.fastq$R1zip[nb_seq], list.fastq$R1filt[nb_seq], list.fastq$R2zip[nb_seq], list.fastq$R2filt[nb_seq],
                     truncQ=3,trimLeft=c(0,0),maxEE = c(10,10),minLen = 100,multithread=T,verbose =F)

print(paste("Done filter&trim",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "filter&trim"

#plot qual before and after
#plotQualityProfile(list.fastq$R1zip[1:4],aggregate =T)
#plotQualityProfile(list.fastq$R1filt[1:4],aggregate =T)

#error rates of sequences
error_R1 <- learnErrors(list.fastq$R1filt[c(1:6,97:102,193:198)],multithread = TRUE)
error_R2 <- learnErrors(list.fastq$R2filt[c(1:6,97:102,193:198)],multithread = TRUE)
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

###get rid of it if less than 2 samples identified
filter = rep(T,length(nb_seq))
for(i in 1:length(nb_seq))
{
  if(length(dada_R2[[i]][1][[1]]) < 2) filter[i] = F
}

dada_R1 = dada_R1[filter]
dada_R2 = dada_R2[filter] 
derep_R1 = derep_R1[filter]
derep_R2 = derep_R2[filter]

sample.names.nb_seq = sample.names[nb_seq]
if(length(filter[filter == F])>0) print(paste("I got rid of sample ",sample.names.nb_seq[filter==F],", It contained a single ASVs...",sep = ""))


#merging of PE
merged <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2,minOverlap = 20, maxMismatch = 3)

print(paste("Done merging",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "merging"

#Identify ASVs
asv <- makeSequenceTable(merged)

print(paste("Done ID ASVs",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "asv"

#Remove chimeras
asv.nochim <- removeBimeraDenovo(asv, method = "consensus", multithread = F,verbose = TRUE) 
asv.nochim.bin <- ifelse(asv.nochim>0,1,0) 

asv.nochim.simplified = asv.nochim
colnames(asv.nochim.simplified) = paste("ASV",c(1:ncol(asv.nochim)),sep = "")

print(paste("Done remove chimeras",Sys.time()))
time=c(time,Sys.time()); names(time)[length(time)] = "chimeras"

#summary stats 
getN <- function(x) sum(getUniques(x))
summary.stats <- cbind(out[filter,], # input & filtered
               round(((out[filter,2]/out[filter,1])*100),2), # % (filtered / input) 
               sapply(merged, getN), # Merged 
               round(((sapply(merged, getN)/out[filter,1])*100),2), # % (Merged / Input)
               rowSums(asv.nochim),# Non-chimeric
               round(((rowSums(asv.nochim)/out[filter,1])*100),2),# % (Non-chimeric / Input)
               rowSums(asv.nochim.bin)) # Number of ASVs per sample 

colnames(summary.stats) <- c("Input","Filter","%Filt/In","Merge","%Mer/In","Non-chim","%Non-chim/In","#ASV/sample") # Column names
rownames(summary.stats) <- sample.names[nb_seq] # Row names
summary.stats<-as.data.frame(summary.stats)
print(summary.stats)

time=c(time,Sys.time()); names(time)[length(time)] = "summary_stats"

#assign taxonomy
asv.taxo <- assignTaxonomy(asv.nochim,refFasta = reference.taxonomy,minBoot = 80, multithread=TRUE)
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
write.table(asv.nochim.simplified,paste("results/asv.",group,sep = ""))
write.table(asv.taxo,paste("results/asv.taxo.",group,sep = ""))



