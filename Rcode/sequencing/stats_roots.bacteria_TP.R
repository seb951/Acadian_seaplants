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
asv = read.table("results/asv/asv.bacteria_root",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F, na.strings = "noNAs")
taxo = read.table("results/asv/asv.taxo.bacteria_root",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F, na.strings = "noNAs")

###REMOVE MITO AND BACTERIAL ASVs...
cp_vector=rep(T,nrow(taxo))
cp_vector[taxo[,4] == "Chloroplast"] = FALSE
cp_vector[taxo[,5] == "Mitochondria"] = FALSE
print("what fraction of reads are mito/cp?")
1- sum(asv[,cp_vector==T]) / sum(asv)

asv = asv[,cp_vector==T]
taxo = taxo[cp_vector==T,]


#tidy up ----
#remove the bad samples (below 4 standard deviations...).
asv.sum = rowSums(asv)
min(asv.sum);max(asv.sum)
asv.sum.mean = mean(asv.sum) 
asv.sum.sd = sd(asv.sum)
#asv.filt = asv[asv.sum > (asv.sum.mean - 4*asv.sum.sd),]
asv.filt = asv[asv.sum > 1000,]
dim(asv.filt)


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

#tomato and pepper specific analyses
design.keep.pepper = design.keep[design.keep[,4]== "pepper",]
design.keep.tomato = design.keep[design.keep[,4]== "tomato",]

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

#dev.new()
#par(mar=c(16,4,4,2))
#barplot(asv.filt.abundants.norm.barplot,beside = F,font = 3, axisnames = TRUE,ylab = "Relative ASV abundance",col = c(cols25(),rep("grey",650)), las = 3, xpd = TRUE)
#legend(0.1,-0.22,fill = cols25()[1:20],legend =  barplot.taxo[1:20],cex = 0.75,xpd =TRUE)
#dev.print(device=pdf, "figures/bacteria/Figure4br_ASVabundance.pdf", onefile=FALSE)
#dev.off()

###barplot by family ----
#summarize with dplyr by family
asv.filt.abundants.norm.FAMILY = as.data.frame(asv.filt.abundants.norm.barplot) %>% group_by(taxo.abundants$Family) %>% summarise_all(sum)

write.table(asv.filt.abundants.norm.FAMILY,"results/asv/asv.filt.abundants.norm.FAMILY_br")


###alpha diversity ----
#prepare a matrix with alpha diversity as "invsimpson" index
asv.filt.abundants.norm.alpha = cbind((vegan::diversity(asv.filt.abundants.norm, index= "invsimpson"))^(1/2),design.keep)
colnames(asv.filt.abundants.norm.alpha)[1] = "alpha"

#linear mixed effect model on alpha diversity (block & replicate are random, replicate is nested in bloc)
lmm.alpha.interactions.tomato <- lme(alpha~fertilization,data = asv.filt.abundants.norm.alpha[design.keep[,4]== "tomato",],random = ~1|bloc/replicate, method = "REML")
lmm.alpha.interactions.pepper <- lme(alpha~fertilization,data = asv.filt.abundants.norm.alpha[design.keep[,4]== "pepper",],random = ~1|bloc/replicate, method = "REML")


anova(lmm.alpha.interactions.tomato)
anova(lmm.alpha.interactions.pepper)
#                       numDF denDF   F-value p-value
#(Intercept)               1    54 729.5858  <.0001
#fertilization             1    54  10.0908  0.0025
#species                   1    54   4.4958  0.0386
#fertilization:species     1    54   0.3433  0.5604

shapiro.test(lmm.alpha.interactions.tomato$residuals[,1]) #normaly distributed SQUARE ROOT
shapiro.test(lmm.alpha.interactions.pepper$residuals[,1]) 

#boxplot it
dev.new(width=7, height=7,units = "cm",noRStudioGD = TRUE)
par(mfrow = c(1,1),mar = c(6,6,4,2))
par(mgp = c(3.5,2,1),cex.lab = 1.5)
names = c("control","amendment","control","amendment")

#
x=boxplot(asv.filt.abundants.norm.alpha$alpha[design.keep[,4]== "tomato"]~asv.filt.abundants.norm.alpha$fertilization[design.keep[,4]== "tomato"],ylab = expression(italic(alpha)~~diversity ("Inverse Simpson Index")),
          main = "Roots - Bacteria",xpd = T,font = 2,xlim = c(0.5,4.5),names = c("",""),ylim = c(2,8))
text(labels = c("a","b"),x = c(1.2,2.2),y = x$stats[4,]+max(x$stats[4,])/30,cex = 1.5,font =3)
mtext(text = c("Tomato","Pepper"),at = c(1.5,3.5),side = 1,line = -23,cex = 1.5)

x=boxplot(asv.filt.abundants.norm.alpha$alpha[design.keep[,4]== "pepper"]~asv.filt.abundants.norm.alpha$fertilization[design.keep[,4]== "pepper"],add = T,at = 3:4,names = c("",""))
#text(labels = c("a","b"),x = c(3.2,4.2),y = x$stats[4,]+max(x$stats[5,])/30,cex = 1.5,font =3)

abline(v = 2.5,lty = 2)
mtext(names,side = 1,at = c(1:4),cex =0.8,line = 1.5,font=3)


dev.print(device=pdf, "figures/bacteria/Figure4br_alpha.pdf", onefile=FALSE)
dev.off()

###PERMANOVA: are community different (beta diversity) according to sampling design  ----
#note that results are essentially the same irrespective of prior normalization (OTU.norm) or not (OTU).
asv.filt.abundants.norm.species = asv.filt.abundants.norm[design.keep[,4]== "tomato",]

#standardization using hellinger transform
asv.filt.abundants.norm.hel <-decostand(asv.filt.abundants.norm.species, "hel")

#PERMANOVA - TOMATO
#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
permanova = adonis(formula=asv.filt.abundants.norm.hel~fertilization,strata=(design.keep.tomato$bloc/design.keep.tomato$replicate), data=design.keep.tomato, permutations=9999, method="bray")
permanova$aov.tab$comparison = "root_bacteria_tomato"
permanova
write.table(permanova$aov.tab,"results/asv/permanova.root_bacteria_tomato")

###PERMANOVA: are community different (beta diversity) according to sampling design  ----
#note that results are essentially the same irrespective of prior normalization (OTU.norm) or not (OTU).
asv.filt.abundants.norm.species = asv.filt.abundants.norm[design.keep[,4]== "pepper",]

#standardization using hellinger transform
asv.filt.abundants.norm.hel <-decostand(asv.filt.abundants.norm.species, "hel")

#PERMANOVA - TOMATO
#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
permanova = adonis(formula=asv.filt.abundants.norm.hel~fertilization,strata=(design.keep.pepper$bloc/design.keep.pepper$replicate), data=design.keep.pepper, permutations=9999, method="bray")
permanova
permanova$aov.tab$comparison = "root_bacteria_pepper"
write.table(permanova$aov.tab,"results/asv/permanova.root_bacteria_pepper")
