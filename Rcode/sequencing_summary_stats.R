
setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants/results") 

library(dplyr)

#getting summary data & clean it up
summary.stats.fungi_root = read.table("asv/summary.stats.fungi_root",header = T,sep = "",stringsAsFactors = F)
summary.stats.bacteria_root = read.table("asv/summary.stats.bacteria_root",header = T,sep = "",stringsAsFactors = F)
summary.stats.bacteria_soil = read.table("asv/summary.stats.bacteria_soil",header = T,sep = "",stringsAsFactors = F)
summary.stats.fungi_soil = read.table("asv/summary.stats.fungi_soil",header = T,sep = "",stringsAsFactors = F)

#summary stat it
summary.stats = rbind(summary.stats.bacteria_root,summary.stats.bacteria_soil,summary.stats.fungi_root,summary.stats.fungi_soil)
summary.stats$treatment =  c(rep("bacteria-root",95),rep("bacteria-soil",192),
                        rep("fungi-root",96),rep("fungi-soil",192))

#
summary.table = summary.stats %>% group_by(treatment) %>% summarise(Nb_samples = n(),
Nb_seq_sumX10e3 = round(sum(Input)/10000),Nb_seq_mean = round(mean(Input)),Nb_seq_mean_filtered = round(mean(Filter)),
Nb_seq_mean_filt_merged = round(mean(Merge)),Nb_seq_mean_filt_merg_non_chimeras = round(mean(Non.chim)),
ASV_persample = 0)

summary.table$ASV_sum = 0 
summary.table$Nb_samples_trimmed = 0
summary.table$ASV_sum_trimmed = 0

for(t in 1:4)
{
asv = read.table(paste("asv/asv.",gsub("-","_",unique(summary.stats$treatment)[t]),sep = ""),row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)
taxo = read.table(paste("asv/asv.taxo.",gsub("-","_",unique(summary.stats$treatment)[t]),sep = ""),row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)

#remove the bad samples (below 4 standard deviations...).
asv.sum = rowSums(asv)
min(asv.sum);max(asv.sum)
asv.sum.mean = mean(asv.sum) 
asv.sum.sd = sd(asv.sum)
#asv.filt = asv[asv.sum > (asv.sum.mean - 4*asv.sum.sd),]
asv.filt = asv[asv.sum > 1000,]

summary.table$ASV_sum[t] = ncol(asv.filt)
summary.table$Nb_samples_trimmed[t] = nrow(asv.filt)


if(i >3) {
  ###REMOVE MITO AND BACTERIAL ASVs...
  cp_vector=rep(T,ncol(asv))
  cp_vector[taxo[,4] == "Chloroplast"] = FALSE
  cp_vector[taxo[,5] == "Mitochondria"] = FALSE
  print("what fraction of reads are mito/cp?")
  print(1- sum(asv.filt[,cp_vector==T]) / sum(asv))
  asv.filt = asv.filt[,cp_vector==T]
}

asv.singletons = rep(0,ncol(asv.filt)) #in how many samples is an ASV found...
for(i in 1:ncol(asv.filt))
{
  asv.singletons[i] = length(asv.filt[asv.filt[,i]!=0,i])
}
asv.filt.abundants = asv.filt[,asv.singletons>(c(5,10,5,10)[t])]
summary.table$ASV_sum_trimmed[t] = ncol(asv.filt.abundants)
print(sum(asv.filt.abundants)/sum(asv.filt))


#No ASV trimmed _persample (mean) 
ASV_persample = rep(0,nrow(asv.filt.abundants))
for(i in 1:nrow(asv.filt)){ASV_persample[i] = length(asv.filt.abundants[i,asv.filt.abundants[i,]>0])}
summary.table$ASV_persample[t] = round(mean(ASV_persample),0)
}

write.table(summary.table, "summary.table",row.names = F)


###table 2 - permanova --- 
permanova.soil_fungi_tomato = read.table("asv/permanova.soil_fungi_tomato",header = T)
permanova.root_fungi_tomato = read.table("asv/permanova.root_fungi_tomato",header = T)
permanova.soil_bacteria_tomato = read.table("asv/permanova.soil_bacteria_tomato",header = T)
permanova.root_bacteria_tomato = read.table("asv/permanova.root_bacteria_tomato",header = T)
permanova.soil_fungi_pepper = read.table("asv/permanova.soil_fungi_pepper",header = T)
permanova.root_fungi_pepper = read.table("asv/permanova.root_fungi_pepper",header = T)
permanova.soil_bacteria_pepper = read.table("asv/permanova.soil_bacteria_pepper",header = T)
permanova.root_bacteria_pepper = read.table("asv/permanova.root_bacteria_pepper",header = T)

permanova.summary = matrix(NA,nrow = 3,ncol = 8)
colnames(permanova.summary) = c("soil_fungi_tomato","soil_fungi_pepper","root_fungi_tomato","root_fungi_pepper","soil_bact_tomato","soil_bact_pepper","root_bact_tomato","root_bact_pepper")
rownames(permanova.summary) = rownames(permanova.soil_fungi_tomato)[1:3]

permanova.summary[,1] = paste(round(permanova.soil_fungi_tomato$R2[1:3],2)," (",round(permanova.soil_fungi_tomato$Pr..F.[1:3],4),")",sep = "")
permanova.summary[,2] = paste(round(permanova.soil_fungi_pepper$R2[1:3],2)," (",round(permanova.soil_fungi_pepper$Pr..F.[1:3],4),")",sep = "")

permanova.summary[1,3] = paste(round(permanova.root_fungi_tomato$R2[1],2)," (",round(permanova.root_fungi_tomato$Pr..F.[1],4),")",sep = "")
permanova.summary[1,4] = paste(round(permanova.root_fungi_pepper$R2[1],2)," (",round(permanova.root_fungi_pepper$Pr..F.[1],4),")",sep = "")

permanova.summary[,5] = paste(round(permanova.soil_bacteria_tomato$R2[1:3],2)," (",round(permanova.soil_bacteria_tomato$Pr..F.[1:3],4),")",sep = "")
permanova.summary[,6] = paste(round(permanova.soil_bacteria_pepper$R2[1:3],2)," (",round(permanova.soil_bacteria_pepper$Pr..F.[1:3],4),")",sep = "")

permanova.summary[1,7] = paste(round(permanova.root_bacteria_tomato$R2[1],2)," (",round(permanova.root_bacteria_tomato$Pr..F.[1],4),")",sep = "")
permanova.summary[1,8] = paste(round(permanova.root_bacteria_pepper$R2[1],2)," (",round(permanova.root_bacteria_pepper$Pr..F.[1],4),")",sep = "")

write.table(permanova.summary,"permanova.summary_PT")






