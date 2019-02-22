

setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants")

#combining figures
library(RColorBrewer)
library(vegan)
library(multipanelfigure)
library(magrittr)
library(dplyr)
#library(DECIPHER) #this is a bioconductor packages
#library(phangorn)
library(phyloseq)

#Figure 4 - barplots -----
for(group in c("bacteria","fungi"))
{
  abr = ifelse(group=="bacteria","b","f")
  #get data and merge it
  asv.filt.abundants.norm.FAMILY_r = read.table(paste("results/asv/asv.filt.abundants.norm.FAMILY_",abr,"r",sep=""),stringsAsFactors = F,na.strings = "none")
  asv.filt.abundants.norm.FAMILY_s = read.table(paste("results/asv/asv.filt.abundants.norm.FAMILY_",abr,"s",sep=""),stringsAsFactors = F,na.strings = "none")

if(group=="bacteria"){
#combine NA + unidentified
asv.filt.abundants.norm.FAMILY_s[asv.filt.abundants.norm.FAMILY_s[,1]=='Unknown_Family',-1] =
asv.filt.abundants.norm.FAMILY_s[asv.filt.abundants.norm.FAMILY_s[,1]=='Unknown_Family',-1] + asv.filt.abundants.norm.FAMILY_s[asv.filt.abundants.norm.FAMILY_s[,1]=='NA',-1]
asv.filt.abundants.norm.FAMILY_s = asv.filt.abundants.norm.FAMILY_s[asv.filt.abundants.norm.FAMILY_s[,1]!='NA',]

asv.filt.abundants.norm.FAMILY_r[asv.filt.abundants.norm.FAMILY_r[,1]=='Unknown_Family',-1] =
asv.filt.abundants.norm.FAMILY_r[asv.filt.abundants.norm.FAMILY_r[,1]=='Unknown_Family',-1] + asv.filt.abundants.norm.FAMILY_r[asv.filt.abundants.norm.FAMILY_r[,1]=='NA',-1]
asv.filt.abundants.norm.FAMILY_r = asv.filt.abundants.norm.FAMILY_r[asv.filt.abundants.norm.FAMILY_r[,1]!='NA',]
}

if(group=="fungi"){
  #combine NA + unidentified
  asv.filt.abundants.norm.FAMILY_r[asv.filt.abundants.norm.FAMILY_r[,1]=='unidentified',-1] =
  asv.filt.abundants.norm.FAMILY_r[asv.filt.abundants.norm.FAMILY_r[,1]=='unidentified',-1] + asv.filt.abundants.norm.FAMILY_r[asv.filt.abundants.norm.FAMILY_r[,1]=='NA',-1]
  asv.filt.abundants.norm.FAMILY_r = asv.filt.abundants.norm.FAMILY_r[asv.filt.abundants.norm.FAMILY_r[,1]!='NA',]
  
  asv.filt.abundants.norm.FAMILY_s[asv.filt.abundants.norm.FAMILY_s[,1]=='unidentified',-1] =
  asv.filt.abundants.nsorm.FAMILY_s[asv.filt.abundants.norm.FAMILY_s[,1]=='unidentified',-1] + asv.filt.abundants.norm.FAMILY_s[asv.filt.abundants.norm.FAMILY_s[,1]=='NA',-1]
  asv.filt.abundants.norm.FAMILY_s = asv.filt.abundants.norm.FAMILY_s[asv.filt.abundants.norm.FAMILY_s[,1]!='NA',]
}

#merge
asv.filt.abundants.norm.FAMILY = merge(asv.filt.abundants.norm.FAMILY_s,asv.filt.abundants.norm.FAMILY_r,by.x = "taxo.abundants.Family", by.y = "taxo.abundants.Family",all =T)

#prepare ggplot object
temp = asv.filt.abundants.norm.FAMILY[order(rowSums(asv.filt.abundants.norm.FAMILY[,2:11],na.rm=T),decreasing = T),]
temp[is.na(temp[,1]),1] = "unknown"
asv.filt.abundants.norm.FAMILY.top10 = temp[1:20,]
barplot.data = as.data.frame(unlist(asv.filt.abundants.norm.FAMILY.top10[,2:11]))
colnames(barplot.data)[1] = "fraction"
n = nrow(asv.filt.abundants.norm.FAMILY.top10)
barplot.data$factor = factor(c(rep("control",n),rep("amendment",n),rep("tomato",n),rep("pepper",n),rep("non-planted",n),rep("planted",n),rep("amendment",n),rep("control",n),rep("tomato",n),rep("pepper",n)),levels=c('control','amendment','tomato','pepper','non-planted',"planted"))
barplot.data$taxonomy = unlist(rep(asv.filt.abundants.norm.FAMILY.top10[,1],ncol(asv.filt.abundants.norm.FAMILY.top10)-1))
barplot.data$soilroot = factor(c(rep("soil",20*6), rep("root",20*4)),levels = c('soil','root'))

###ggplot
dev.new()
x = colorRampPalette(brewer.pal(12,"Paired"))
p=ggplot() + labs(title = group,fill = "Taxonomy (family)",y ="Relative abundance of ASVs", x = "") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_bar(aes(y = fraction, x = factor, fill = taxonomy),
           data = barplot.data,stat="identity") +
  scale_fill_manual(values = x(nrow(asv.filt.abundants.norm.FAMILY.top10))) +
 facet_grid(rows=vars(soilroot))
p 

dev.print(device=pdf, paste("figures/Figure4_FAMILY_barplots_",group,".pdf",sep=""), onefile=FALSE)
dev.off()
}

#Figure 5 - alpha -----
figure4 = multi_panel_figure(width = 14,height = 14, unit = "inch",rows = 2,columns = 2)

figure4 %<>% fill_panel("figures/fungi/Figure4fs_alpha.pdf",row = 1, column = 1)
figure4 %<>% fill_panel("figures/fungi/Figure4fr_alpha.pdf",row = 1, column = 2)
figure4 %<>% fill_panel("figures/bacteria/Figure4bs_alpha.pdf",row = 2, column = 1)
figure4 %<>% fill_panel("figures/bacteria/Figure4br_alpha.pdf",row = 2, column = 2)
figure4 %>% save_multi_panel_figure(filename = "figures/Figure4_alpha.pdf")






