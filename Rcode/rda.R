#CCAs

###building the canonical correspondence analyses + picking out candidate genes

library(vegan)
library(dplyr)
setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants")


#Figure5 - CCA fungi and Figure6 - CCA bacteria-----
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
#fungi, then bacteria
for(spp in c("f","b")) 
{
#soil vs root loop
for(biotope in c("soil","root"))
{
  if(biotope == "soil") asv.filt.abundants.norm = read.table(paste("results/asv.filt.abundants.norm_",spp,"s",sep=""))
  if(biotope == "soil") taxo.abundants = read.table(paste("results/taxo.abundants_",spp,"s",sep=""))
  if(biotope == "root") asv.filt.abundants.norm = read.table(paste("results/asv.filt.abundants.norm_",spp,"r",sep=""))
  if(biotope == "root") taxo.abundants = read.table(paste("results/taxo.abundants_",spp,"r",sep=""))
  
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
    rda = rda(asv.filt.abundants.norm.species.hel, productivity.norm.keep.species[,c(6,8,9,10)])
    
    #verify model significance
    print(anova.cca(rda))
    
    #plot
    if(species == "Tomato" & biotope == "soil") dev.new()
    if(species == "Tomato" & biotope == "soil") dev.new(width=7, height=7,units = "inch",noRStudioGD=T)
    if(species == "Tomato" & biotope == "soil") par(mar=c(5,5,2,1),xpd =T,mgp = c(3.5,2,1),mfrow=c(2,2))
    
    rda.plot = plot(rda,scaling = "species",type = "n",font=2,cex.lab=1.2,font.lab=2,main = paste(species," (",biotope," - ",ifelse(spp == "b","bacteria","fungi"),")",sep = ""))
    points(rda,scaling = "species",display = c("sp"),pch = 3,col = "red",cex = 0.6)
    points(rda,scaling = "species",display = c("bp"),pch = 3,col = "blue")
    text(rda,scaling = "species",display = c("bp"),col = "blue",cex = 0.8,font = 2)
    
    col = ifelse(productivity.norm.keep.species[,2] == "F+","goldenrod4","darkgrey")
    
    #add plot sites
    text(rda.plot$sites,labels = rownames(rda.plot$sites),cex = 0.7, col = col,font = 2,adj = 0.8)
    
    #Candidate ASVs (top10?) closest to arrowheads (excluding avg fruit weigth)
    factors = c(1,3,4)     #remove avg. fruit weight. it is orthogonal to the other variables
    arrow_x = mean(rda.plot$biplot[factors,1]*ordiArrowMul(rda,display = "bp"))
    arrow_y = mean(rda.plot$biplot[factors,2]*ordiArrowMul(rda,display = "bp"))
    
    dist_x = arrow_x - rda.plot$species[,1]
    dist_y = arrow_y - rda.plot$species[,2]
    dist = abs(dist_x)+abs(dist_y)
    
    #candidate top10 and plot it.
    candidate.top10 = rda.plot$species[order(dist),][1:10,]
    points(candidate.top10,pch = 3, lwd =1,col = "blue")
    
    candidate.ASV = rbind(candidate.ASV,cbind(candidate.top10,species))
    
    
    #get abundance data for the candidate ten most associated with productivity
    out=NULL
    for(i in 1:10){out = cbind(out,asv.filt.abundants.norm.species[,colnames(asv.filt.abundants.norm.species) == rownames(candidate.top10)[i]])}
    
#    as.data.frame(out)
    
    #variance partinioning (4 productivity measures in addition to fertilisation treatment)
    varpart = varpart(
      asv.filt.abundants.norm.species,~.,
      ifelse(productivity.norm.keep.species[,c(2)]=="F+",1,0),
      data= productivity.norm.keep.species[,c(6,8,9,10)],transfo="hel")
    #plot
    #plot(varpart,alpha = 50,bg = c("orange","blue"),Xnames = c("Fertilisation","Productivity"))
    
  }
  
  #candidates with taxonomy
  candidate.ASV = cbind(rownames(candidate.ASV),candidate.ASV)
  colnames(candidate.ASV)[1] = "ASV"
  
  #taxo with ASV
  taxo.abundants$ASV = rownames(taxo.abundants)
  
  #inner join
  candidate.ASV.taxo = inner_join(as.data.frame(candidate.ASV,stringsAsFactors = F),taxo.abundants,by = "ASV")
  
  #save candidates
  if(biotope == "soil") write.table(candidate.ASV.taxo,paste("results/candidate.ASV.",spp,"s.txt",sep=""),row.names = T, col.names = T, quote = T)
  if(biotope == "root") write.table(candidate.ASV.taxo,paste("results/candidate.ASV.",spp,"r.txt",sep=""),row.names = T, col.names = T, quote = T)
  
  if(biotope == "root" & spp == "f") dev.print(device=pdf,"figures/fungi/Figure5f_RDA.pdf", onefile=FALSE)
  if(biotope == "root" & spp == "b") dev.print(device=pdf,"figures/bacteria/Figure6b_RDA.pdf", onefile=FALSE)
  if(biotope == "root") dev.off()
    }
  }

#dev.print(device=pdf,"figures/varpart.pdf", onefile=FALSE)
#dev.off()                                                                            
                                                                                             