setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants")

###setting things up ----
library(nlme)
library(vegan)
library(pals)
library(ape)
library(dplyr)
library(ggplot2)
#load ASVS and taxonomy information ----
dev.new()
asv = read.table("results/asv/asv.fungi_root",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)
taxo = read.table("results/asv/asv.taxo.fungi_root",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)

#tidy up ----
#remove the bad samples (below 4 standard deviations...).
asv.sum = rowSums(asv)
min(asv.sum);max(asv.sum)
asv.sum.mean = mean(asv.sum) 
asv.sum.sd = sd(asv.sum)
asv.filt = asv[asv.sum > (asv.sum.mean - 4*asv.sum.sd),]
asv.filt = asv[asv.sum > 1000,]


#98% of ASVs are present in less than 5 individuals (5% of samples): this is huge. and justifies the clustering I think for soil bacteria....
asv.singletons = rep(0,ncol(asv.filt))
for(i in 1:ncol(asv.filt))
{
  asv.singletons[i] = length(asv.filt[asv.filt[,i]!=0,i])
}
asv.filt.abundants = asv.filt[,asv.singletons>5]
taxo.abundants = taxo[asv.singletons>5,]

#what fraction of ASV do you got rid off
print(length(asv.singletons[asv.singletons<5])/length(asv.singletons))  

#what percentage of reads did you keep
sum(asv.filt.abundants)/sum(asv.filt)


#relative abundance
asv.filt.abundants.norm = asv.filt.abundants/rowSums(asv.filt.abundants)

###statistal analyses on ASVs----
design = read.table("reference_material/experimental_design.txt", header = TRUE, stringsAsFactors = FALSE)

#keep samples of interest in a new design.keep dataframe.
design.keep = NULL
for(i in 1:nrow(asv.filt.abundants.norm))
{
  design.keep = rbind(design.keep,design[design[,2] == rownames(asv.filt.abundants.norm)[i],])
}

asv.filt.abundants.norm.design = cbind(asv.filt.abundants.norm, design.keep)

#select and ASV (here ASV1)
ASV1 = asv.filt.abundants.norm.design[,c(1,(ncol(asv.filt.abundants.norm.design)-9):ncol(asv.filt.abundants.norm.design))]

#linear mixed effect model  (block & replicate are random, replicate is nested in bloc)/ for ASV1 ----
lmm1 <- lme(ASV1~fertilization+species,data = ASV1,random = ~1|bloc/replicate, method = "REML")
anova(lmm1)

#               numDF denDF  F-value p-value
#(Intercept)       1    55 596.9096  <.0001
#fertilization     1    55   4.8193  0.0324
#species           1    55  12.3121  0.0009

###barplots of all ASVs ----
###to look at diversity among fertilization, soil/roots, planted, species----
#first group things according to treatment and species.
asv.filt.abundants.norm.barplot = cbind(colMeans(asv.filt.abundants.norm[design.keep$fertilization == 0, ]), 
                         colMeans(asv.filt.abundants.norm[design.keep$fertilization == 1, ]),
                         colMeans(asv.filt.abundants.norm[design.keep$species == "tomato", ]), 
                         colMeans(asv.filt.abundants.norm[design.keep$species == "pepper", ]))

colnames(asv.filt.abundants.norm.barplot) = c("fertilization", "control","tomato","pepper")      

#barplot (TOP20)
barplot.taxo = c(1:nrow(taxo))
for(i in 1:length(barplot.taxo))
{
  barplot.taxo[i]=paste(taxo[i,2:6],collapse = ";")
}

dev.new()
par(mar=c(16,4,4,2))
barplot(asv.filt.abundants.norm.barplot,beside = F,font = 3, axisnames = TRUE,ylab = "Relative ASV abundance",col = c(cols25(),rep("grey",10000)), las = 3, xpd = TRUE)
legend(0.1,-0.22,fill = cols25()[1:20],legend =  barplot.taxo[1:20],cex = 0.75,xpd =TRUE)
dev.print(device=pdf, "figures/fungi/Figure4fr_ASVabundance.pdf", onefile=FALSE)
dev.off()

#summarize with dplyr.
#summarize with dplyr by FAMILY
asv.filt.abundants.norm.FAMILY = as.data.frame(asv.filt.abundants.norm.barplot) %>% group_by(taxo.abundants$Family) %>% summarise_all(sum)
write.table(asv.filt.abundants.norm.FAMILY,"results/asv/asv.filt.abundants.norm.FAMILY_fr")

###alpha diversity ----
#prepare a matrix with alpha diversity as "invsimpson" index
asv.filt.abundants.norm.alpha = cbind((vegan::diversity(asv.filt.abundants.norm, index= "invsimpson"))^(1/2),design.keep)
colnames(asv.filt.abundants.norm.alpha)[1] = "alpha"

#linear mixed effect model on alpha diversity (block & replicate are random, replicate is nested in bloc)
lmm.alpha.interactions <- lme(alpha~fertilization*species,data = asv.filt.abundants.norm.alpha,random = ~1|bloc/replicate, method = "REML")
anova(lmm.alpha.interactions)

#                       numDF denDF   F-value p-value
#(Intercept)               1    54 729.5858  <.0001
#fertilization             1    54  10.0908  0.0025
#species                   1    54   4.4958  0.0386
#fertilization:species     1    54   0.3433  0.5604

shapiro.test(lmm.alpha.interactions$residuals[,1]) #normaly distributed SQUARE ROOT
#W = 0.9427, p-value = 0.001243

#boxplot it?
#boxplot it
dev.new(width=7, height=7,units = "cm",noRStudioGD = TRUE)
par(mfrow = c(1,1),mar = c(5,4,4,2))
names = expression(italic(control),italic(fertilized),italic(Pepper),italic(Tomato))

#
x=boxplot(asv.filt.abundants.norm.alpha$alpha~asv.filt.abundants.norm.alpha$fertilization,ylab = expression(italic(alpha)~~diversity ("Inverse Simpson Index")),
        main = "Roots - Fungi",xpd = T,font = 2,xlim = c(0.5,4.5),names = c("",""))
text(labels = c("a","b"),x = c(1.2,2.2),y = x$stats[4,]+max(x$stats[4,])/30,cex = 1.5,font =3)

x=boxplot(asv.filt.abundants.norm.alpha$alpha~asv.filt.abundants.norm.alpha$species,add = T,at = 3:4,names = c("",""))
text(labels = c("a","b"),x = c(3.2,4.2),y = x$stats[4,]+max(x$stats[5,])/30,cex = 1.5,font =3)

abline(v = 2.5,lty = 2)
mtext(names,side = 1,at = c(1:4),cex =0.8,line = 1)

dev.print(device=pdf, "figures/fungi/Figure4fr_alpha.pdf", onefile=FALSE)
dev.off()

###PcoA ----
###PERMANOVA: are community different (beta diversity) according to sampling design  ----
#note that results are essentially the same irrespective of prior normalization (OTU.norm) or not (OTU).

#standardization using hellinger transform
asv.filt.abundants.norm.hel <-decostand(asv.filt.abundants.norm, "hel")

#PERMANOVA
#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
permanova = adonis(formula=asv.filt.abundants.norm.hel~fertilization*species,strata=(design.keep$bloc/design.keep$replicate), data=design.keep, permutations=9999, method="bray")
permanova$aov.tab$comparison = "root_fungi"
write.table(permanova$aov.tab,"results/asv/permanova.root_fungi")


#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#fertilization          1    0.8643 0.86429  10.345 0.08261 0.0001           
#species                1    2.7476 2.74764  32.886 0.26264 0.0001           
#fertilization:species  1    0.4165 0.41646   4.985 0.03981 0.0015           
#Residuals             77    6.4333 0.08355         0.61494                  
#Total                 80   10.4617                 1.00000          
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

####PcoA: good to illustrate the community structure----
#Calculating Bray-Curtis dissimilarity matrix on the hellinger transformed data
asv.filt.abundants.norm.hel.bray <-vegdist(asv.filt.abundants.norm.hel, method="bray")

#Calculating PCoA
asv.filt.abundants.norm.hel.bray.pcoa<-pcoa(asv.filt.abundants.norm.hel.bray)

#How many axes represent more variability (17)
bs = asv.filt.abundants.norm.hel.bray.pcoa$values$Broken_stick
length(bs[bs>mean(bs)])

#PVE of first 2 axes (4.7% & 3.8%)
axis.1.2 = round(c(asv.filt.abundants.norm.hel.bray.pcoa$values$Relative_eig[1],asv.filt.abundants.norm.hel.bray.pcoa$values$Relative_eig[2]),4)*100

#Ploting the PCoAs - with fertilization as empty circles
#crops are "darkred","darkblue","darkorange
col = design.keep$species
col = gsub("tomato","darkorange",col);col=gsub("pepper","darkblue",col)
pch.pcoa = rep(0,nrow(design.keep))
pch.pcoa[design.keep$fertilization == 1] = 19; pch.pcoa[design.keep$fertilization == 0] = 21

dev.new()
plot(asv.filt.abundants.norm.hel.bray.pcoa$vectors[,1],asv.filt.abundants.norm.hel.bray.pcoa$vectors[,2],col = col, pch = pch.pcoa,
     ylab = paste("PC2",sep = ""), xlab = paste("PC1",sep = ""))
abline(h = 0,lty=3,cex =0.7);abline(v=0,lty=3,cex =0.7)
legend(0.3,0.2,fill = c("darkorange","darkblue"),legend = c("  Tomato","  Pepper"),box.lwd = 1)
legend(0.36,0.2,fill = rep("transparent",3), border = c("darkorange","darkblue"),legend = rep("",2),box.lwd = 0,box.col = "transparent")

#inoculation text (this is a pain...)
text(c(0.42,0.49),c(0.22,0.22),c("fertilized","control"),srt = 45,pos = 3,font =3)
dev.print(device=pdf, "figures/fungi/Figure5fr_pcoa.pdf", onefile=FALSE)
dev.off()



###RDA-----
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
#####

#keep samples of interest in a new design.keep dataframe.
productivity.norm.keep = NULL
for(i in 1:nrow(asv.filt.abundants.norm))
{
  productivity.norm.keep = rbind(productivity.norm.keep,productivity.norm[productivity.norm[,15] == rownames(asv.filt.abundants.norm)[i],])
}

candidate.ASV = NULL
#Do an RDA for TOMATOES or PEPPERS only (given that we know that both crops differ a lot...)
for(species in c("Tomato","Pepper"))
{
  #keep only one species in the ASV matrix and the PRODUCTIVITY matrix
  productivity.norm.keep.species = productivity.norm.keep[productivity.norm.keep[,1] == species,]
  asv.filt.abundants.norm.species = asv.filt.abundants.norm[productivity.norm.keep[,1] == species,]
  
  #keep only the non-zeros...
  asv.filt.abundants.norm.species = asv.filt.abundants.norm.species[,colSums(asv.filt.abundants.norm.species)!=0]
  
  #hellinger transform
  asv.filt.abundants.norm.species.hel = decostand(asv.filt.abundants.norm.species, "hel")
  
  #RDA (Constrained Ordination)
  rda = cca(asv.filt.abundants.norm.species.hel, productivity.norm.keep.species[,c(6,8,9,10)])
  
  #verify model significance
  print(anova.cca(rda))
  
  #plot
  dev.new(width=7, height=7,units = "inch")
  par(mar=c(6,6,4,2),xpd =T,mgp = c(3.5,2,1))
  rda.plot = plot(rda,scaling = "species",type = "n")
  xlims = c(min(rda.plot$sites[,1])*1.1,max(rda.plot$sites[,1])*1.1)
  ylims = c(min(rda.plot$sites[,2])*1.1,max(rda.plot$sites[,2])*1.1)
  rda.plot2=plot(rda,scaling = "species",display = c("sp", "bp"),main = paste(species," (root - fungi)",sep = ""),xlim = xlims,ylim = ylims)
  #  xlim = c(-4,4),ylim = c(-4,4)
  #coloring for fertilization
  col = ifelse(productivity.norm.keep.species[,2] == "F+","goldenrod4","darkgrey")
  
  #add plot sites
  text(rda.plot$sites,labels = rownames(rda.plot$sites),cex = 0.7, col = col,font = 2,adj = 0.8)
  
  #Candidate ASVs (top10?) closest to arrowheads (excluding avg fruit weigth)
  factors = c(1,3,4)     #remove avg. fruit weight for tomato
  arrow_x = mean(rda.plot2$biplot[factors,1]*attr(rda.plot2$biplot,"arrow.mul"))
  arrow_y = mean(rda.plot2$biplot[factors,2]*attr(rda.plot2$biplot,"arrow.mul"))
  
  dist_x = arrow_x - rda.plot2$species[,1]
  dist_y = arrow_y - rda.plot2$species[,2]
  dist = abs(dist_x)+abs(dist_y)
  
  #candidate top10 and plot it.
  candidate.top10 = rda.plot$species[order(dist),][1:10,]
  points(candidate.top10,lwd =4,col = "red")
  
  candidate.ASV = rbind(candidate.ASV,cbind(candidate.top10,species))
  
  dev.print(device=pdf, paste("figures/fungi/Figure6fr_RDA_",species,".pdf",sep = ""), onefile=FALSE)
  dev.off()
}

#candidates with taxonomy
candidate.ASV = cbind(rownames(candidate.ASV),candidate.ASV)
colnames(candidate.ASV)[1] = "ASV"

#taxo with ASV
taxo.abundants$ASV = rownames(taxo.abundants)

#inner join
candidate.ASV.taxo = inner_join(as.data.frame(candidate.ASV,stringsAsFactors = F),taxo.abundants,by = "ASV")

#
write.table(candidate.ASV.taxo,"results/candidate.ASV.fr.txt",row.names = T, col.names = T, quote = T)




###sandbox -----
#this is essentially to show that the aproach above (taking the top10 candidates) makes sense. As you get the same pattern, with a cor.test approach on the productivity variables....
if(1==2){
cortest = data.frame(colnames(asv.filt.abundants.norm.species.hel) ,stringsAsFactors = F)
colnames(cortest)[1] = "asv"
cortest$pvalue = 0
cortest$cor = 0
for( i in 1:ncol(asv.filt.abundants.norm.species.hel))
{
  cortest$pvalue[i] = cor.test(asv.filt.abundants.norm.species.hel[,i],productivity.norm.keep.species[,10])$p.value
  cortest$cor[i] = cor.test(asv.filt.abundants.norm.species.hel[,i],productivity.norm.keep.species[,10])$estimate
}
head(cortest[order(cortest$cor,decreasing = T),],10)
}