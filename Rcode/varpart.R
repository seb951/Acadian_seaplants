

####variance partionning

####but how many species to use????
###show the first 15 (they represent 50% of the species abundance!)
###what proportion of variance in the productivity is due to fertilization versus presence of ASV? 

library(vegan)
library(dplyr)
setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants")






#Figure5 - CCA fungi -----
#fungi
####setting things up
productivity = read.table("results/plant_productivity_data.tsv", header = T, stringsAsFactors = F,sep = "\t")
productivity.norm = productivity

#we do a sqrt transformation on all variables to help in normalizing (it helps in all cases, but it's not perfect. Some residuals are still not normally distributed. But the effect is so strong that we don't care too much.)
productivity.norm$fruit.number = sqrt(productivity$fruit.number)
productivity.norm$fruits.weight = sqrt(productivity$fruits.weight)
productivity.norm$shoots.fresh.weight = sqrt(productivity$shoots.fresh.weight)
productivity.norm$shoots.dry.weight = sqrt(productivity$shoots.dry.weight)
productivity.norm$roots.fresh.weight = sqrt(productivity$roots.fresh.weight)
productivity.norm$roots.dry.weight = sqrt(productivity$roots.dry.weight)
#####de

#soil vs root loop
for(biotope in c("soil","root"))
{
  if(biotope == "soil") asv.filt.abundants.norm = read.table("results/asv.filt.abundants.norm_fs")
  if(biotope == "soil") taxo.abundants = read.table("results/taxo.abundants_fs")
  if(biotope == "root") asv.filt.abundants.norm = read.table("results/asv.filt.abundants.norm_fr")
  if(biotope == "root") taxo.abundants = read.table("results/taxo.abundants_fr")
  
  #keep samples of interest in a new design.keep dataframe.
  if(biotope == "soil") 
  {
    productivity.norm.keep = productivity.norm[0,]
    for(i in 1:nrow(asv.filt.abundants.norm))
    {
      temp = productivity.norm[productivity.norm[,16] == rownames(asv.filt.abundants.norm)[i],]
      if(nrow(temp)==1) productivity.norm.keep[i,] = temp
      if(nrow(temp)==0) productivity.norm.keep[i,] = rep(0,16)
    }
  }
  
  if(biotope == "root") 
  {  
    #keep samples of interest in a new design.keep dataframe.
    productivity.norm.keep = NULL
    for(i in 1:nrow(asv.filt.abundants.norm))
    {
      productivity.norm.keep = rbind(productivity.norm.keep,productivity.norm[productivity.norm[,15] == rownames(asv.filt.abundants.norm)[i],])
    }
  }
  
  candidate.ASV = NULL
  #Do an RDA for TOMATOES or PEPPERS only (given that we know that both crops differ a lot...)
  for(species in c("Tomato","Pepper"))
  {
    #keep only one species in the ASV matrix and the PRODUCTIVITY matrix
    productivity.norm.keep.species = productivity.norm.keep[productivity.norm.keep[,1] == species,]
    asv.filt.abundants.norm.species = asv.filt.abundants.norm[productivity.norm.keep[,1] == species,]
    
    #hellinger transform
    asv.filt.abundants.norm.species.hel = decostand(asv.filt.abundants.norm.species, "hel")
    
    #RDA (Constrained Ordination)
    rda = cca(asv.filt.abundants.norm.species.hel, productivity.norm.keep.species[,c(6,8,9,10)])
    
    
    
    
#get abundance data for the candidate ten most associated with productivity
out=NULL
for(i in 1:10){out = cbind(out,asv.filt.abundants.norm.species[,colnames(asv.filt.abundants.norm.species) == rownames(candidate.top10)[i]])}

#variance partinioning (4 productivity measures in addition to fertilisation treatment)
varpart = varpart(
  as.data.frame(out),~.,
  ifelse(productivity.norm.keep.species[,c(2)]=="F+",1,0),
  data= productivity.norm.keep.species[,c(6,8,9,10)],transfo="hel")
#plot
plot(varpart,alpha = 50,bg = c("orange","blue"))

