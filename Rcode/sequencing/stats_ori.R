#!/usr/local/bin/Rscript

setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants")

###setting things up ----
library(nlme)
library(vegan)
library(pals)
library(ape)

#load ASVS and taxonomy information ----
asv = read.table("results/asv.bacteria_soil",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)
taxo = read.table("results/asv.taxo.bacteria_soil",row.names = 1, sep = " ",header =TRUE,stringsAsFactors = F)


#tidy up ----
#remove the bad samples (below 4 standard deviations...).
asv.sum = rowSums(asv)
min(asv.sum);max(asv.sum)
asv.sum.mean = mean(asv.sum) 
asv.sum.sd = sd(asv.sum)
asv.filt = asv[asv.sum > (asv.sum.mean - 4*asv.sum.sd),]

#98% of ASVs are present in less than 10 individuals (5% of samples): this is huge. and justifies the clustering I think for soil bacteria....
asv.singletons = rep(0,ncol(asv.filt))
for(i in 1:ncol(asv.filt))
{
  asv.singletons[i] = length(asv.filt[asv.filt[,i]!=0,i])
}
print(length(asv.singletons[asv.singletons<10])/length(asv.singletons))  

asv.filt[,asv.singletons<10])


#Remove ASVs that represent less than 1% of the sequences
asv.cumul = rep(0,ncol(asv.filt)+1)
for(i in 1:ncol(asv.filt))
  {
    asv.percent = sum((colSums(asv.filt)/sum(asv.filt))[i])
    asv.cumul[i+1] = asv.cumul[i]+asv.percent
  }
#95% or 99% of the sequences...
asv.filt.99 = asv.filt[,asv.cumul<0.99]
asv.filt.95 = asv.filt[,asv.cumul<0.95]
nrow(asv.filt)
nrow(asv.filt.99)
nrow(asv.filt.95)

#relative abundance
asv.filt.95.norm = asv.filt.95/rowSums(asv.filt.95)

###statistal analyses on ASVs----
design = read.table("reference_material/experimental_design.txt", header = TRUE, stringsAsFactors = FALSE)

#keep samples of interest in a new design.keep dataframe.
design.keep = NULL
for(i in 1:nrow(asv.filt.95.norm))
{
  design.keep = rbind(design.keep,design[design[,2] == rownames(asv.filt.95.norm)[i],])
}

asv.filt.95.norm.design = cbind(asv.filt.95.norm, design.keep)

#select and ASV (here ASV1)
ASV1 = asv.filt.95.norm.design[,c(1,(ncol(asv.filt.95.norm.design)-9):ncol(asv.filt.95.norm.design))]

#linear mixed effect model (block is random) for ASV1 ----
lmm1 <- lme(ASV1~fertilization+replicate,data = ASV1,random = ~1|bloc, method = "REML")
anova(lmm1)

###barplots of all ASVs ----
###to look at diversity among fertilization, soil/roots, planted, species----
#first group things according to treatment and species.
asv.filt.95.norm.barplot = cbind(colMeans(asv.filt.95.norm[design.keep$fertilization == 0, ]), 
                         colMeans(asv.filt.95.norm[design.keep$fertilization == 1, ]),
                         colMeans(asv.filt.95.norm[design.keep$species == "tomato", ]), 
                         colMeans(asv.filt.95.norm[design.keep$species == "pepper", ]),
                         colMeans(asv.filt.95.norm[design.keep$planted == 0, ]), 
                         colMeans(asv.filt.95.norm[design.keep$planted == 1, ]),
                         colMeans(asv.filt.95.norm[design.keep$type == "soil", ]), 
                         colMeans(asv.filt.95.norm[design.keep$type == "root", ]))

colnames(asv.filt.95.norm.barplot) = c("fertilization", "no.fertiliz.","tomato","pepper","non-planted","planted","soil","roots")      

#barplot (TOP10)
barplot.taxo = c(1:nrow(taxo))
for(i in 1:length(barplot.taxo))
{
  barplot.taxo[i]=paste(taxo[i,2:6],collapse = ";")
}

dev.new()
par(mar=c(16,4,4,2))
barplot(asv.filt.95.norm.barplot[1:10,],beside = F,font = 3, axisnames = TRUE,ylab = "Relative ASV abundance",col = c(cols25(),rep("grey",60)), las = 3, xpd = TRUE)
legend(0.1,-0.22,fill = cols25()[1:10],legend =  barplot.taxo[1:10],cex = 0.75,xpd =TRUE)
dev.print(device=pdf, "figures/figure1_OTUabundance.pdf", onefile=FALSE)
dev.off()


###alpha diversity ----
#prepare a matrix with alpha diversity as "invsimpson" index
asv.filt.95.norm.alpha = cbind(sqrt(diversity(asv.filt.95.norm, index= "invsimpson")),design.keep)
colnames(asv.filt.95.norm.alpha)[1] = "alpha"

#linear mixed effect model on alpha diversity (block is random)
lmm.alpha <- lme(alpha~fertilization+replicate,data = asv.filt.95.norm.alpha,random = ~1|bloc, method = "REML")
anova(lmm.alpha)

shapiro.test(lmm.alpha$residuals) #normaly distributed (otherwise can do a log or sqrt?)

#boxplot it?
dev.new(width=14, height=6,units = "cm",noRStudioGD = TRUE)
par(mfrow = c(1,4),mar = c(6,5,4,2))

#fertilization
boxplot(asv.filt.95.norm.alpha$alpha~asv.filt.95.norm.alpha$fertilization,beside = F,font = 3, axisnames = T,ylab = expression(italic(alpha)~~diversity),las = 3, xpd = T,main = "Fertilization",cex.main = 1.5,cex.lab = 1.5,names = expression(italic(control),italic(fertilized)))

#species
boxplot(asv.filt.95.norm.alpha$alpha~asv.filt.95.norm.alpha$fertilization,beside = F,font = 3, axisnames = T,ylab = expression(italic(alpha)~~diversity),las = 3, xpd = T,main = "Species",cex.main = 1.5,cex.lab = 1.5,names = expression(italic(tomato),italic(pepper)))

#soil/roots
boxplot(asv.filt.95.norm.alpha$alpha~asv.filt.95.norm.alpha$fertilization,beside = F,font = 3, axisnames = T,ylab = expression(italic(alpha)~~diversity),las = 3, xpd = T,main = "Soils/roots",cex.main = 1.5,cex.lab = 1.5,names = expression(italic(soil),italic(roots)))

#planted? 
boxplot(asv.filt.95.norm.alpha$alpha~asv.filt.95.norm.alpha$fertilization,beside = F,font = 3, axisnames = T,ylab = expression(italic(alpha)~~diversity),las = 3, xpd = T,main = "Planted?",cex.main = 1.5,cex.lab = 1.5,names = expression(italic(control),italic(fertilized)))

dev.print(device=pdf, "figures/figure2_alpha.pdf", onefile=FALSE)
dev.off()

###PcoA ----
###PERMANOVA: are community different (beta diversity) according to sampling design  ----
#note that results are essentially the same irrespective of prior normalization (OTU.norm) or not (OTU).

#standardization using hellinger transform
asv.filt.95.norm.hel <-decostand(asv.filt.95.norm, "hel")

#PERMANOVA (I only keep the species*growing_stage interaction as the other interactions were all NS)
#The only assumption of PERMANOVA is independence of samples (I think, but could be wrong here)
adonis(formula=asv.filt.95.norm.hel ~ fertilization+replicate+bloc, data=design.keep, permutations=9999, method="bray")

#               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#fertilization  1   0.18142 0.181419  4.3171 0.08704 0.0001 ***
#replicate      1   0.04873 0.048727  1.1595 0.02338 0.2336    
#bloc           1   0.04728 0.047283  1.1252 0.02268 0.2667    
#Residuals     43   1.80701 0.042023         0.86690           
#Total         46   2.08444                  1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#### posthoc permanova to look at species + growing stage effect----
#beta diversity per species to look at growth stage?----

#use species AND growing stage as factors (9999 perms for more precision given that one p is close to 0.05) 
#pairwise.perm.manova(dist(OTU,method = "euclidean"),fact = paste(design.keep[,1],design.keep[,2],sep = "_"),nperm = 9999)

####PcoA: good to illustrate the community structure----
#Calculating Bray-Curtis dissimilarity matrix on the hellinger transformed data
asv.filt.95.norm.hel.bray <-vegdist(asv.filt.95.norm.hel, method="bray")

#Calculating PCoA
asv.filt.95.norm.hel.bray.pcoa<-pcoa(dist(asv.filt.95.norm.hel.bray))

#How many axes represent more variability (17)
bs = asv.filt.95.norm.hel.bray.pcoa$values$Broken_stick
length(bs[bs>mean(bs)])

#PVE of first 2 axes (4.7% & 3.8%)
axis.1.2 = round((asv.filt.95.norm.hel.bray.pcoa$values$Broken_stick/sum(asv.filt.95.norm.hel.bray.pcoa$values$Broken_stick))[1:2],4)*100

#Ploting the PCoAs - with inoculation as empty circles
#crops are "darkred","darkblue","darkorange
col = design.keep$species
col = gsub("tomato","darkorange",col);col=gsub("pepper","darkblue",col)
pch.pcoa = rep(0,nrow(design.keep))
pch.pcoa[design.keep$fertilization == 1] = 19; pch.pcoa[design.keep$fertilization == 0] = 21

dev.new()
plot(asv.filt.95.norm.hel.bray.pcoa$vectors[,1],asv.filt.95.norm.hel.bray.pcoa$vectors[,2],col = col, pch = pch.pcoa,
     ylab = paste("PC2",sep = ""), xlab = paste("PC1",sep = ""))
abline(h = 0,lty=3,cex =0.7);abline(v=0,lty=3,cex =0.7)
legend(0.3,0.2,fill = c("darkorange","darkblue"),legend = c("  Tomato","  Pepper"),box.lwd = 1)
legend(0.36,0.2,fill = rep("transparent",3), border = c("darkorange","darkblue"),legend = rep("",2),box.lwd = 0,box.col = "transparent")

#inoculation text (this is a pain...)
text(c(0.42,0.49),c(0.22,0.22),c("fertilized","control"),srt = 45,pos = 3,font =3)
dev.print(device=pdf, "figures/figure3_pcoa.pdf", onefile=FALSE)
dev.off()


###sandbox -----
dev.new()
par(mfrow = c(4,5))
for(i in 1:20)
{
  ax1 = sample(c(97:158),1)
    ax2 = sample(c(97:158),1)
  plot.default(asv[ax1,],asv[ax2,],cex.lab = 0.5,xlab = ax1,ylab = ax2)
  }

