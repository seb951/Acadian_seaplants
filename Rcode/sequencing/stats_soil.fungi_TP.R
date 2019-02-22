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
asv = read.table("results/asv/asv.fungi_soil",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)
taxo = read.table("results/asv/asv.taxo.fungi_soil",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)

#tidy up ----
#remove the bad samples (below 4 standard deviations... OR 1000 reads...).
asv.sum = rowSums(asv)
min(asv.sum);max(asv.sum)
asv.sum.mean = mean(asv.sum) 
asv.sum.sd = sd(asv.sum)
asv.filt = asv[asv.sum > (asv.sum.mean - 4*asv.sum.sd),]
asv.filt = asv[asv.sum > 1000,]


#98% of ASVs are present in less than 10 individuals (5% of samples): this is huge. and justifies the clustering I think for soil bacteria....
asv.singletons = rep(0,ncol(asv.filt))
for(i in 1:ncol(asv.filt))
{
  asv.singletons[i] = length(asv.filt[asv.filt[,i]!=0,i])
}
asv.filt.abundants = asv.filt[,asv.singletons>10]
taxo.abundants = taxo[asv.singletons>10,]


#what fraction of ASV do you got rid off
print(length(asv.singletons[asv.singletons<10])/length(asv.singletons))  

#what percentage of reads did you keep
sum(asv.filt.abundants)/sum(asv.filt)

#relative abundance
asv.filt.abundants.norm = asv.filt.abundants/rowSums(asv.filt.abundants)

#save it 
write.table(asv.filt.abundants.norm,"results/asv.filt.abundants.norm_fs")
write.table(taxo.abundants,"results/taxo.abundants_fs")

###statistal analyses on ASVs----
design = read.table("reference_material/experimental_design.txt", header = TRUE, stringsAsFactors = FALSE)

#keep samples of interest in a new design.keep dataframe.
design.keep = NULL
for(i in 1:nrow(asv.filt.abundants.norm))
{
  design.keep = rbind(design.keep,design[design[,2] == rownames(asv.filt.abundants.norm)[i],])
}

#tomato and pepper specific analyses
design.keep.pepper = design.keep[design.keep[,4]== "pepper",]
design.keep.tomato = design.keep[design.keep[,4]== "tomato",]


asv.filt.abundants.norm.design = cbind(asv.filt.abundants.norm, design.keep)

#select and ASV (here ASV1)
ASV1 = asv.filt.abundants.norm.design[,c(1,(ncol(asv.filt.abundants.norm.design)-9):ncol(asv.filt.abundants.norm.design))]

#linear mixed effect model  (block & replicate are random, replicate is nested in bloc)/ for ASV1 ----
lmm1 <- lme(ASV1~fertilization+species+planted,data = ASV1,random = ~1|bloc/replicate, method = "REML")
anova(lmm1)
#numDF denDF   F-value p-value
#(Intercept)       1   162 112.98950  <.0001
#fertilization     1   162   0.49873  0.4811
#species           1   162   0.02763  0.8682
#planted           1   162   2.53759  0.1131

###barplots of all ASVs ----
###to look at diversity among fertilization, soil/roots, planted, species----
#first group things according to treatment and species.
asv.filt.abundants.norm.barplot = cbind(colMeans(asv.filt.abundants.norm[design.keep$fertilization == 0, ]), 
                                        colMeans(asv.filt.abundants.norm[design.keep$fertilization == 1, ]),
                                        colMeans(asv.filt.abundants.norm[design.keep$species == "tomato", ]), 
                                        colMeans(asv.filt.abundants.norm[design.keep$species == "pepper", ]),
                                        colMeans(asv.filt.abundants.norm[design.keep$planted == 0, ]), 
                                        colMeans(asv.filt.abundants.norm[design.keep$planted == 1, ]))

colnames(asv.filt.abundants.norm.barplot) = c("control", "fertilization","tomato","pepper","control","planted")      

#barplot (TOP10)
barplot.taxo = c(1:nrow(taxo.abundants))
for(i in 1:length(barplot.taxo))
{
  barplot.taxo[i]=paste(taxo.abundants[i,2:6],collapse = ";")
}

#dev.new()
#par(mar=c(16,4,4,2))
#barplot(asv.filt.abundants.norm.barplot,beside = F,font = 3, axisnames = TRUE,cex.names = 0.6, cex.lab = 0.6,ylab = "Relative ASV abundance",col = c(cols25(),rep("grey",10000)), las = 3, xpd = TRUE)
#legend(0.1,-0.12,fill = cols25()[1:20],legend =  barplot.taxo[1:20],cex = 0.75,xpd =TRUE)
#dev.print(device=pdf, "figures/fungi/Figure4fs_ASVabundance.pdf", onefile=FALSE)
#dev.off()

#summarize with dplyr.
asv.filt.abundants.norm.FAMILY = as.data.frame(asv.filt.abundants.norm.barplot) %>% group_by(taxo.abundants$Family) %>% summarise_all(sum)
barplot.data = as.data.frame(asv.filt.abundants.norm.FAMILY,row.names = taxo.abundants$Family)
barplot.data = barplot.data[order(rowSums(barplot.data[,2:5]),decreasing = T),]

dev.new()
par(mar=c(16,4,4,2))
barplot(as.matrix(barplot.data)[,2:7],beside = F,font = 3, axisnames = TRUE,ylab = "Relative ASV abundance",col = c(cols25(),rep("grey",650)), las = 3, xpd = TRUE)
legend(0.1,-0.32,fill = cols25(),legend =  barplot.data[,1],cex = 0.75,xpd =TRUE)


#summarize with dplyr.
#summarize with dplyr by FAMILY
asv.filt.abundants.norm.FAMILY = as.data.frame(asv.filt.abundants.norm.barplot) %>% group_by(taxo.abundants$Family) %>% summarise_all(sum)
#write.table(asv.filt.abundants.norm.FAMILY,"results/asv/asv.filt.abundants.norm.FAMILY_fs")


###alpha diversity ----
#prepare a matrix with alpha diversity as "invsimpson" index
asv.filt.abundants.norm.alpha = cbind((vegan::diversity(asv.filt.abundants.norm, index= "invsimpson"))^(1/1),design.keep)
colnames(asv.filt.abundants.norm.alpha)[1] = "alpha"

#linear mixed effect model on alpha diversity (block & replicate are random, replicate is nested in bloc)
lmm.alpha.interactions.tomato <- lme(alpha~fertilization*planted,data = asv.filt.abundants.norm.alpha[design.keep[,4]== "tomato",],random = ~1|bloc/replicate, method = "REML")
lmm.alpha.interactions.pepper <- lme(alpha~fertilization*planted,data = asv.filt.abundants.norm.alpha[design.keep[,4]== "pepper",],random = ~1|bloc/replicate, method = "REML")

anova(lmm.alpha.interactions.tomato)
anova(lmm.alpha.interactions.pepper)
#                       numDF denDF   F-value p-value
#(Intercept)               1    54 729.5858  <.0001
#fertilization             1    54  10.0908  0.0025
#species                   1    54   4.4958  0.0386
#fertilization:species     1    54   0.3433  0.5604

shapiro.test(lmm.alpha.interactions.tomato$residuals[,1]) #normaly distributed SQUARE ROOT
shapiro.test(lmm.alpha.interactions.pepper$residuals[,1]) 

####Note, here I tested a full model comprisingly both the root and soil fungi. This is to show that there is more diversity in soil than roots...
####asv.filt.abundants.norm.alpha.root.soil are a rbind of the root and soil fungi.
####lmm.alpha.interactions.pepper <- lme(alpha~fertilization*type*species,data = asv.filt.abundants.norm.alpha.root.soil,random = ~1|bloc/replicate, method = "REML")
####anova(lmm.alpha.interactions.pepper)
#numDF denDF   F-value p-value
#(Intercept)                    1   239 1073.2288  <.0001
#fertilization                  1   239    0.0065  0.9360
#type                           1   239  899.4884  <.0001
#species                        1   239   12.3809  0.0005
#fertilization:type             1   239    0.0437  0.8345
#fertilization:species          1   239    2.9713  0.0860
#type:species                   1   239    4.9069  0.0277
#fertilization:type:species     1   239    1.4220  0.2343

#boxplot it
dev.new(width=7, height=7,units = "cm",noRStudioGD = TRUE)
par(mfrow = c(1,1),mar = c(6,6,4,2))
par(mgp = c(3.5,2,1),cex.lab = 1.5)
names = c("control","amendment","control","planted","control","amendment","control","planted")

#
x=boxplot(asv.filt.abundants.norm.alpha$alpha[design.keep[,4]== "tomato"]~asv.filt.abundants.norm.alpha$fertilization[design.keep[,4]== "tomato"],ylab = expression(italic(alpha)~~diversity ("Inverse Simpson Index")),
          main = "Soil - Fungi",xpd = T,font = 2,xlim = c(0.5,8.5),names = c("",""),ylim = c(0,50))
#text(labels = c("a","b"),x = c(1.2,2.2),y = x$stats[4,]+max(x$stats[4,])/30,cex = 1.5,font =3)
mtext(text = c("Tomato","Pepper"),at = c(2.5,6.5),side = 1,line = -23,cex = 1.5)

x=boxplot(asv.filt.abundants.norm.alpha$alpha[design.keep[,4]== "tomato"]~asv.filt.abundants.norm.alpha$planted[design.keep[,4]== "tomato"],add = T,at = 3:4,names = c("",""))
text(labels = c("a","b"),x = c(3.2,4.2),y = x$stats[4,]+max(x$stats[5,])/30,cex = 1.5,font =3)

x=boxplot(asv.filt.abundants.norm.alpha$alpha[design.keep[,4]== "pepper"]~asv.filt.abundants.norm.alpha$fertilization[design.keep[,4]== "pepper"],add = T,at = 5:6,names = c("",""))
#text(labels = c("a","b"),x = c(5.2,6.2),y = x$stats[4,]+max(x$stats[5,])/30,cex = 1.5,font =3)

x=boxplot(asv.filt.abundants.norm.alpha$alpha[design.keep[,4]== "pepper"]~asv.filt.abundants.norm.alpha$planted[design.keep[,4]== "pepper"],add = T,at = 7:8,names = c("",""))
text(labels = c("a","b"),x = c(7.2,8.2),y = x$stats[4,]+max(x$stats[5,])/30,cex = 1.5,font =3)
abline(v = 4.5,lty = 2)
mtext(names,side = 1,at = c(1:8),cex =0.8,line = 1.4,font=3)

dev.print(device=pdf, "figures/fungi/Figure4fs_alpha.pdf", onefile=FALSE)
dev.off()

###PERMANOVA: are community different (beta diversity) according to sampling design  ----
#note that results are essentially the same irrespective of prior normalization (OTU.norm) or not (OTU).
asv.filt.abundants.norm.species = asv.filt.abundants.norm[design.keep[,4]== "tomato",]

#standardization using hellinger transform
asv.filt.abundants.norm.hel <-decostand(asv.filt.abundants.norm.species, "hel")

#PERMANOVA - TOMATO
#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
permanova = adonis(formula=asv.filt.abundants.norm.hel~fertilization*planted,strata=(design.keep.tomato$bloc/design.keep.tomato$replicate), data=design.keep.tomato, permutations=9999, method="bray")
permanova$aov.tab$comparison = "soil_fungi_tomato"
permanova
write.table(permanova$aov.tab,"results/asv/permanova.soil_fungi_tomato")

###PERMANOVA: are community different (beta diversity) according to sampling design  ----
#note that results are essentially the same irrespective of prior normalization (OTU.norm) or not (OTU).
asv.filt.abundants.norm.species = asv.filt.abundants.norm[design.keep[,4]== "pepper",]

#standardization using hellinger transform
asv.filt.abundants.norm.hel <-decostand(asv.filt.abundants.norm.species, "hel")

#PERMANOVA - TOMATO
#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
permanova = adonis(formula=asv.filt.abundants.norm.hel~fertilization*planted,strata=(design.keep.pepper$bloc/design.keep.pepper$replicate), data=design.keep.pepper, permutations=9999, method="bray")
permanova
permanova$aov.tab$comparison = "soil_fungi_pepper"
write.table(permanova$aov.tab,"results/asv/permanova.soil_fungi_pepper")


