
setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants")

###setting things up ----
library(nlme)
library(vegan)
library(pals)
library(ape)
library(ggplot2)
library(dplyr)

#load ASVS and taxonomy information ----
dev.new()
asv = read.table("results/asv/asv.bacteria_soil",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)
taxo = read.table("results/asv/asv.taxo.bacteria_soil",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)

#tidy up ----
#remove the bad samples (below 4 standard deviations...).
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

#linear mixed effect model (block is random) for ASV1 ----
lmm1 <- lme(ASV1~fertilization+species+planted,data = ASV1,random = ~1|bloc/replicate, method = "REML")
anova(lmm1)

#             numDF denDF   F-value p-value
#(Intercept)       1   165 203.32592  <.0001
#fertilization     1   165   2.87183  0.0920
#species           1   165   2.93756  0.0884
#planted           1   165   1.85800  0.1747

###barplots of all ASVs ----
###to look at diversity among fertilization, soil/roots, planted, species----
#first group things according to treatment and species.
asv.filt.abundants.norm.barplot = cbind(colMeans(asv.filt.abundants.norm[design.keep$fertilization == 0, ]), 
                         colMeans(asv.filt.abundants.norm[design.keep$fertilization == 1, ]),
                         colMeans(asv.filt.abundants.norm[design.keep$species == "tomato", ]), 
                         colMeans(asv.filt.abundants.norm[design.keep$species == "pepper", ]),
                         colMeans(asv.filt.abundants.norm[design.keep$planted == 0, ]), 
                         colMeans(asv.filt.abundants.norm[design.keep$planted == 1, ]))

colnames(asv.filt.abundants.norm.barplot) = c("control", "fertilization","tomato","pepper","non-planted","planted")      

#barplot
barplot.taxo = c(1:nrow(taxo.abundants))
for(i in 1:length(barplot.taxo))
{
  barplot.taxo[i]=paste(taxo.abundants[i,2:6],collapse = ";")
}

dev.new()
par(mar=c(16,4,4,2))
barplot(asv.filt.abundants.norm.barplot,beside = F,font = 3, axisnames = TRUE,cex.names = 0.6, cex.lab = 0.6,ylab = "Relative ASV abundance",col = c(cols25(),rep("grey",11000)), las = 3, xpd = TRUE)
legend(0.1,-0.12,fill = cols25()[1:20],legend =  barplot.taxo[1:20],cex = 0.75,xpd =TRUE)
dev.print(device=pdf, "figures/bacteria/Figure4bs_ASVabundance.pdf", onefile=FALSE)
dev.off()

#summarize with dplyr.
#summarize with dplyr by FAMILY
asv.filt.abundants.norm.FAMILY = as.data.frame(asv.filt.abundants.norm.barplot) %>% group_by(taxo.abundants$Family) %>% summarise_all(sum)

#top10
temp = asv.filt.abundants.norm.FAMILY[order(rowSums(asv.filt.abundants.norm.FAMILY[,2:7]),decreasing = T),]
temp[is.na(temp),1] = "unknown"
asv.filt.abundants.norm.FAMILY.top10 = temp[1:10,]
barplot.data = as.data.frame(unlist(asv.filt.abundants.norm.FAMILY.top10[,2:7]))
colnames(barplot.data)[1] = "fraction"
n = nrow(asv.filt.abundants.norm.FAMILY.top10)
barplot.data$treatment =c(rep("control",n),rep("fertilization",n),rep("tomato",n),rep("pepper",n),rep("non-planted",n),rep("planted",n))
barplot.data$taxonomy = unlist(rep(asv.filt.abundants.norm.FAMILY.top10[,1],ncol(asv.filt.abundants.norm.FAMILY.top10)-1))
###ggplot ----
dev.new()
p = ggplot() + labs(title = "Soil - bacteria",fill = "Taxonomy (family)",y ="Relative abundance of ASVs") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_bar(aes(y = fraction, x = treatment, fill = taxonomy),
           data = barplot.data,stat="identity")
p + scale_fill_brewer(palette = "Set3")
dev.print(device=pdf, "figures/bacteria/Figure4bs_FAMILY_ASVabundance.pdf", onefile=FALSE)
dev.off()






###alpha diversity ----
#prepare a matrix with alpha diversity as "invsimpson" index
asv.filt.abundants.norm.alpha = cbind((vegan::diversity(asv.filt.abundants.norm, index= "invsimpson")^(1/1)),design.keep)
colnames(asv.filt.abundants.norm.alpha)[1] = "alpha"

#linear mixed effect model on alpha diversity (bloc/replicate are random, none of the interactions are significant)
lmm.alpha <- lme(alpha~fertilization+species+planted,data = asv.filt.abundants.norm.alpha,random = ~1|bloc/replicate, method = "REML")
anova(lmm.alpha)
#              numDF denDF   F-value p-value
#(Intercept)       1   165 11971.730  <.0001
#fertilization     1   165    17.131  0.0001
#species           1   165     1.891  0.1709
#planted           1   165   138.925  <.0001

shapiro.test(lmm.alpha$residuals[,1]) #fixed effects are  normaly distributed WITHOUT A sqrt
#W = 0.99111, p-value = 0.2834s

#boxplot it
dev.new(width=9, height=7,units = "cm",noRStudioGD = TRUE)
par(mfrow = c(1,1),mar = c(5,5,4,2))
names = expression(italic(control),italic(fertilized),italic(Pepper),italic(Tomato),italic(planted),italic(non-planted))

#
boxplot(asv.filt.abundants.norm.alpha$alpha~asv.filt.abundants.norm.alpha$fertilization,ylab = expression(italic(alpha)~~diversity ("Inverse Simpson Index")),
        main = "Soil - Bacteria",xpd = T,main = "Species",font = 2,xlim = c(0.5,6.5),names = c("",""))
boxplot(asv.filt.abundants.norm.alpha$alpha~asv.filt.abundants.norm.alpha$species,add = T,at = 3:4,names = c("",""))
boxplot(asv.filt.abundants.norm.alpha$alpha~asv.filt.abundants.norm.alpha$planted,add = T,at = 5:6,names = c("",""))

abline(v = c(2.5,4.5),lty = 2)
mtext(names,side = 1,at = c(1:6),cex =0.8,line = 1)
names = expression(italic(control),italic(fertilized),italic(Pepper),italic(Tomato),italic(planted),italic(non-planted))

dev.print(device=pdf, "figures/bacteria/Figure4bs_alpha.pdf", onefile=FALSE)
dev.off()

###PcoA ----
###PERMANOVA: are community different (beta diversity) according to sampling design  ----
#note that results are essentially the same irrespective of prior normalization (OTU.norm) or not (OTU).

#standardization using hellinger transform
asv.filt.abundants.norm.hel <-decostand(asv.filt.abundants.norm, "hel")

#PERMANOVA (I only keep the species*growing_stage interaction as the other interactions were all NS)
#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
permanova = adonis(formula=asv.filt.abundants.norm.hel ~ fertilization*planted*species,strata=(design.keep$bloc/design.keep$replicate), data=design.keep, permutations=9999, method="bray")
permanova$aov.tab$comparison = "soil_bacteria"
write.table(permanova$aov.tab,"results/asv/permanova.soil_bacteria")



#                               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F) comparison
#fertilization                   1    0.9382  0.9382  8.9626 0.03667 0.0001           
#planted                         1    3.2819  3.2819 31.3526 0.12829 0.0001           
#species                         1    0.4054  0.4054  3.8727 0.01585 0.0003           
#fertilization:planted           1    0.5476  0.5476  5.2314 0.02141 0.0001           
#fertilization:species           1    0.6592  0.6592  6.2975 0.02577 0.0001           
#planted:species                 1    0.2882  0.2882  2.7536 0.01127 0.0041           
#fertilization:planted:species   1    0.2005  0.2005  1.9149 0.00784 0.0355           
#Residuals                     184   19.2607  0.1047         0.75291                  
#Total                         191   25.5817                 1.00000         
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#### posthoc permanova to look at species + fertil effects----
#pairwise.perm.manova(dist(OTU,method = "euclidean"),fact = paste(design.keep[,1],design.keep[,2],sep = "_"),nperm = 9999)

####PcoA: good to illustrate the community structure----
#Calculating Bray-Curtis dissimilarity matrix on the hellinger transformed data
asv.filt.abundants.norm.hel.bray <-vegdist(asv.filt.abundants.norm.hel, method="bray")

#Calculating PCoA
asv.filt.abundants.norm.hel.bray.pcoa<-pcoa(dist(asv.filt.abundants.norm.hel.bray))

#How many axes represent more variability (17)
bs = asv.filt.abundants.norm.hel.bray.pcoa$values$Broken_stick
length(bs[bs>mean(bs)])

#PVE of first 2 axes (4.7% & 3.8%)
axis.1.2 = round((asv.filt.abundants.norm.hel.bray.pcoa$values$Broken_stick/sum(asv.filt.abundants.norm.hel.bray.pcoa$values$Broken_stick))[1:2],4)*100

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
dev.print(device=pdf, "figures/bacteria/Figure5bs_pcoa.pdf", onefile=FALSE)
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
productivity.norm.keep = productivity.norm[0,]
for(i in 1:nrow(asv.filt.abundants.norm))
{
  temp = productivity.norm[productivity.norm[,16] == rownames(asv.filt.abundants.norm)[i],]
  if(nrow(temp)==1) productivity.norm.keep[i,] = temp
  if(nrow(temp)==0) productivity.norm.keep[i,] = rep(0,16)
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
  
  #plot
  dev.new(width=7, height=7,units = "inch")
  par(mar=c(6,6,4,2),xpd =T,mgp = c(3.5,2,1))
  rda.plot = plot(rda,scaling = "species",type = "n")
  xlims = c(min(rda.plot$sites[,1])*1.1,max(rda.plot$sites[,1])*1.1)
  ylims = c(min(rda.plot$sites[,2])*1.1,max(rda.plot$sites[,2])*1.1)
  rda.plot2=plot(rda,scaling = "species",display = c("sp", "bp"),main = paste(species," (soil - bacteria)",sep = ""),xlim = xlims,ylim = ylims)
  #  xlim = c(-4,4),ylim = c(-4,4)
  #coloring for fertilization
  col = ifelse(productivity.norm.keep.species[,2] == "F+","goldenrod4","darkgrey")
  
  #add plot sites
  text(rda.plot$sites,labels = rownames(rda.plot$sites),cex = 0.7, col = col,font = 2,adj = 0.8)
  
  #Candidate ASVs (top10?) closest to arrowheads (excluding avg fruit weigth)
  factors = c(1,3,4)     #remove avg. fruit weight
  arrow_x = mean(rda.plot2$biplot[factors,1]*attr(rda.plot2$biplot,"arrow.mul"))
  arrow_y = mean(rda.plot2$biplot[factors,2]*attr(rda.plot2$biplot,"arrow.mul"))
  
  dist_x = arrow_x - rda.plot2$species[,1]
  dist_y = arrow_y - rda.plot2$species[,2]
  dist = abs(dist_x)+abs(dist_y)
  
  #candidate top10 and plot it.
  candidate.top10 = rda.plot$species[order(dist),][1:10,]
  points(candidate.top10,lwd =4,col = "red")
  
  candidate.ASV = rbind(candidate.ASV,cbind(candidate.top10,species))

  dev.print(device=pdf, paste("figures/bacteria/Figure6bs_RDA_",species,".pdf",sep = ""), onefile=FALSE)
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
write.table(candidate.ASV.taxo,"results/candidate.ASV.bs.txt",row.names = T, col.names = T, quote = T)


####sandox ----
if(1==2)
  
{
  rda_t = cca(productivity.norm.keep.species[,c(6,7,9,10)],asv.filt.abundants.norm.species.hel)
  dev.new()
  plot(rda_t)
  mod <- varpart(asv.filt.abundants.norm.species.hel, ~ productivity.norm.keep.species[,2], ~ productivity.norm.keep.species[,6],, data = productivity.norm.keep.species)
  plot(mod)
  
}