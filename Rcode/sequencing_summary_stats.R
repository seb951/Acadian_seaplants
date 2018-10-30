


setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants/results") 

library(dplyr)

#getting summary data & clean it up
summary.stats.fungi_root = read.table("asv/summary.stats.fungi_root",header = T,sep = "",stringsAsFactors = F)
summary.stats.bacteria_root = read.table("asv/summary.stats.bacteria_root",header = T,sep = "",stringsAsFactors = F)
summary.stats.bacteria_soil = read.table("asv/summary.stats.bacteria_soil",header = T,sep = "",stringsAsFactors = F)
summary.stats.fungi_soil = read.table("asv/summary.stats.fungi_soil",header = T,sep = "",stringsAsFactors = F)

#summary stat it:
summary.stats = rbind(summary.stats.bacteria_root,summary.stats.bacteria_soil,summary.stats.fungi_root,summary.stats.fungi_soil)
summary.stats$treatment =  c(rep("bacteria_root",95),rep("bacteria_soil",192),
                        rep("fungi_root",96),rep("fungi_soil",192))

#
summary.table = summary.stats %>% group_by(treatment) %>% summarise(Nb_samples = n(),
Nb_seq_sumX10e3 = round(sum(Input)/10000),Nb_seq_mean = round(mean(Input)),Nb_seq_mean_filtered = round(mean(Filter)),
Nb_seq_mean_filt_merged = round(mean(Merge)),Nb_seq_mean_filt_merg_non_chimeras = round(mean(Non.chim)),
ASV_persample = round(mean(X.ASV.sample)))

summary.table$ASV_sum = 0 
summary.table$Nb_samples_trimmed = 0
summary.table$ASV_sum_trimmed = 0

for(t in 1:4)
{
asv = read.table(paste("asv/asv.",unique(summary.stats$treatment)[t],sep = ""),row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)
#remove the bad samples (below 4 standard deviations...).
asv.sum = rowSums(asv)
min(asv.sum);max(asv.sum)
asv.sum.mean = mean(asv.sum) 
asv.sum.sd = sd(asv.sum)
asv.filt = asv[asv.sum > (asv.sum.mean - 4*asv.sum.sd),]
asv.filt = asv[asv.sum > 1000,]
summary.table$ASV_sum[t] = ncol(asv.filt)
summary.table$Nb_samples_trimmed[t] = nrow(asv.filt)

asv.singletons = rep(0,ncol(asv.filt)) #in how many samples is an ASV found...
for(i in 1:ncol(asv.filt))
{
  asv.singletons[i] = length(asv.filt[asv.filt[,i]!=0,i])
}
asv.filt.abundants = asv.filt[,asv.singletons>(c(5,10,5,10)[t])]
summary.table$ASV_sum_trimmed[t] = ncol(asv.filt.abundants)
print(sum(asv.filt.abundants)/sum(asv.filt))
}
write.table(summary.table, "summary.table",row.names = F)


###table 2 - permanova --- 
permanova.soil_fungi = read.table("results/asv/permanova.soil_fungi",header = T)
permanova.root_fungi = read.table("results/asv/permanova.root_fungi",header = T)
permanova.soil_bacteria = read.table("results/asv/permanova.soil_bacteria",header = T)
permanova.root_bacteria = read.table("results/asv/permanova.root_bacteria",header = T)

permanova.summary = matrix(NA,nrow = 7,ncol = 4)
colnames(permanova.summary) = c("soil_fungi","root_fungi","soil_bact","root_bact")
rownames(permanova.summary) = rownames(permanova.soil_fungi)[1:7]

permanova.summary[,1] = paste(round(permanova.soil_fungi$R2[1:7],2)," (",round(permanova.soil_fungi$Pr..F.[1:7],4),")",sep = "")
permanova.summary[c(1,3,5),2] = paste(round(permanova.root_fungi$R2[1:3],2)," (",round(permanova.root_fungi$Pr..F.[1:3],4),")",sep = "")
permanova.summary[,3] = paste(round(permanova.soil_bacteria$R2[1:7],2)," (",round(permanova.soil_bacteria$Pr..F.[1:7],4),")",sep = "")
permanova.summary[c(1,3,5),4] = paste(round(permanova.root_bacteria$R2[1:3],2)," (",round(permanova.root_bacteria$Pr..F.[1:3],4),")",sep = "")

write.table(permanova.summary,"results/permanova.summary")


permanova.summary






