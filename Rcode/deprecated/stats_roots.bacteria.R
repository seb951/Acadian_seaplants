setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants")

###setting things up ----
library(nlme)
library(vegan)
library(pals)
library(ape)
library(dplyr)
library(ggplot2)
source("/Library/Frameworks/R.framework/Versions/3.4/Resources/library/Legendre_numericalecology/anova.2way.R")
#load ASVS and taxonomy information ----
dev.new()
asv = read.table("results/asv/asv.bacteria_root",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)
taxo = read.table("results/asv/asv.taxo.bacteria_root",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)

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

#what fraction of ASV did you get rid off
print(length(asv.singletons[asv.singletons<5])/length(asv.singletons))  

#what percentage of reads did you keep
sum(asv.filt.abundants)/sum(asv.filt)


#relative abundance
asv.filt.abundants.norm = asv.filt.abundants/rowSums(asv.filt.abundants)

#save it 
write.table(asv.filt.abundants.norm,"results/asv.filt.abundants.norm_br")
write.table(taxo.abundants,"results/taxo.abundants_br")

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

#              numDF denDF   F-value p-value
#(Intercept)       1    68 20102.958  <.0001
#fertilization     1    68    21.868  <.0001
#species           1    68   179.835  <.0001
 


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
barplot(asv.filt.abundants.norm.barplot,beside = F,font = 3, axisnames = TRUE,ylab = "Relative ASV abundance",col = c(cols25(),rep("grey",650)), las = 3, xpd = TRUE)
legend(0.1,-0.22,fill = cols25()[1:20],legend =  barplot.taxo[1:20],cex = 0.75,xpd =TRUE)
dev.print(device=pdf, "figures/bacteria/Figure4br_ASVabundance.pdf", onefile=FALSE)
dev.off()

###barplot by family ----
#summarize with dplyr by family
asv.filt.abundants.norm.FAMILY = as.data.frame(asv.filt.abundants.norm.barplot) %>% group_by(taxo.abundants$Family) %>% summarise_all(sum)

asv.filt.abundants.norm.FAMILY[47,1] = "chloroplast*"
asv.filt.abundants.norm.FAMILY[26,1] = "mitochondria*"

write.table(asv.filt.abundants.norm.FAMILY,"results/asv/asv.filt.abundants.norm.FAMILY_br")


###alpha diversity ----
#prepare a matrix with alpha diversity as "invsimpson" index
asv.filt.abundants.norm.alpha = cbind(vegan::diversity(asv.filt.abundants.norm, index= "invsimpson")^(1/1),design.keep)
colnames(asv.filt.abundants.norm.alpha)[1] = "alpha"

#linear mixed effect model on alpha diversity (block & replicate are random, replicate is nested in bloc)
lmm.alpha.interactions <- lme(alpha~fertilization*species,data = asv.filt.abundants.norm.alpha,random = ~1|bloc/replicate, method = "REML")
anova(lmm.alpha.interactions) #interactions

#                     numDF denDF  F-value p-value
#(Intercept)               1    68 481.0934  <.0001
#fertilization             1    68  17.2742   1e-04
#species                   1    68 359.6885  <.0001
#fertilization:species     1    68  17.5512   1e-04

shapiro.test(lmm.alpha.interactions$residuals) #normaly distributed (otherwise can do a log or sqrt?)
#W = 0.91509, p-value = 1.251e-11

#boxplot it?
dev.new(width=7, height=7,units = "cm",noRStudioGD = TRUE)
par(mfrow = c(1,1),mar = c(5,4,4,2))
names = expression(italic(control),italic(fertilized),italic(Pepper),italic(Tomato))

#
x=boxplot(asv.filt.abundants.norm.alpha$alpha~asv.filt.abundants.norm.alpha$fertilization,ylab = expression(italic(alpha)~~diversity ("Inverse Simpson Index")),
        main = "Roots - Bacteria",xpd = T,main = "Species",font = 2,xlim = c(0.5,4.5),names = c("",""))
text(labels = c("a","b"),x = c(1.2,2.2),y = x$stats[4,]+max(x$stats[4,])/20,cex = 1.5,font =3)

x=boxplot(asv.filt.abundants.norm.alpha$alpha~asv.filt.abundants.norm.alpha$species,add = T,at = 3:4,names = c("",""))
text(labels = c("a","b"),x = c(3.2,4.2),y = x$stats[4,]+max(x$stats[5,])/20,cex = 1.5,font =3)

abline(v = 2.5,lty = 2)
mtext(names,side = 1,at = c(1:4),cex =0.8,line = 1)

dev.print(device=pdf, "figures/bacteria/Figure4br_alpha.pdf", onefile=FALSE)
dev.off()

###PcoA ----
###PERMANOVA: are community different (beta diversity) according to sampling design  ----
#note that results are essentially the same irrespective of prior normalization (OTU.norm) or not (OTU).

#standardization using hellinger transform
asv.filt.abundants.norm.hel <-decostand(asv.filt.abundants.norm, "hel")

#PERMANOVA
#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
permanova=adonis(formula=asv.filt.abundants.norm.hel~fertilization*species,strata=(design.keep$bloc/design.keep$replicate), data=design.keep, permutations=9999, method="bray")
permanova$aov.tab$comparison = "root_bacteria"
write.table(permanova$aov.tab,"results/asv/permanova.root_bacteria")

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#fertilization          1    0.8420  0.8420  16.832 0.06760  1e-04           
#species                1    6.4947  6.4947 129.838 0.52147  1e-04           
#fertilization:species  1    0.5659  0.5659  11.313 0.04544  2e-04           
#Residuals             91    4.5520  0.0500         0.36549                  
#Total                 94   12.4546                 1.00000 
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
dev.print(device=pdf, "figures/bacteria/Figure5br_pcoa.pdf", onefile=FALSE)
dev.off()


###RDA ----

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
  rda.plot2=plot(rda,scaling = "species",display = c("sp", "bp"),main = paste(species," (roots - bacteria)",sep = ""),xlim = xlims,ylim = ylims)
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
  
  dev.print(device=pdf, paste("figures/bacteria/Figure6br_RDA_",species,".pdf",sep = ""), onefile=FALSE)
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
write.table(candidate.ASV.taxo,"results/candidate.ASV.br.txt",row.names = T, col.names = T, quote = T)

