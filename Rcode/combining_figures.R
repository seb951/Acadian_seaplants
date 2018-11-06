

setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants")

#combining figures
library(RColorBrewer)
library(multipanelfigure)
library(magrittr)
library(dplyr)
library(DECIPHER) #this is a bioconductor packages
library(phangorn)
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
asv.filt.abundants.norm.FAMILY_s = asv.filt.abundants.norm.FAMILY_s[-nrow(asv.filt.abundants.norm.FAMILY_s),]
}

if(group=="fungi"){
  #combine NA + unidentified
  asv.filt.abundants.norm.FAMILY_r[asv.filt.abundants.norm.FAMILY_r[,1]=='unidentified',-1] =
    asv.filt.abundants.norm.FAMILY_r[asv.filt.abundants.norm.FAMILY_r[,1]=='unidentified',-1] + asv.filt.abundants.norm.FAMILY_r[asv.filt.abundants.norm.FAMILY_r[,1]=='NA',-1]
  asv.filt.abundants.norm.FAMILY_r = asv.filt.abundants.norm.FAMILY_r[-nrow(asv.filt.abundants.norm.FAMILY_r),]
  
  asv.filt.abundants.norm.FAMILY_s[asv.filt.abundants.norm.FAMILY_s[,1]=='unidentified',-1] =
    asv.filt.abundants.norm.FAMILY_s[asv.filt.abundants.norm.FAMILY_s[,1]=='unidentified',-1] + asv.filt.abundants.norm.FAMILY_s[asv.filt.abundants.norm.FAMILY_s[,1]=='NA',-1]
  asv.filt.abundants.norm.FAMILY_s = asv.filt.abundants.norm.FAMILY_s[-nrow(asv.filt.abundants.norm.FAMILY_s),]
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
barplot.data$treatment = factor(c(rep("control",n),rep("fertilization",n),rep("tomato",n),rep("pepper",n),rep("non-planted",n),rep("planted",n),rep("fertilization",n),rep("control",n),rep("tomato",n),rep("pepper",n)),levels=c('control','fertilization','tomato','pepper','planted','non-planted'))
barplot.data$taxonomy = unlist(rep(asv.filt.abundants.norm.FAMILY.top10[,1],ncol(asv.filt.abundants.norm.FAMILY.top10)-1))
barplot.data$soilroot = factor(c(rep("soil",20*6), rep("root",20*4)),levels = c('soil','root'))

###ggplot
dev.new()
x = colorRampPalette(brewer.pal(12,"Paired"))
p=ggplot() + labs(title = group,fill = "Taxonomy (family)",y ="Relative abundance of ASVs") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_bar(aes(y = fraction, x = treatment, fill = taxonomy),
           data = barplot.data,stat="identity") +
  scale_fill_manual(values = x(nrow(asv.filt.abundants.norm.FAMILY.top10))) +
 facet_grid(rows=vars(soilroot))
p 

dev.print(device=pdf, paste("figures/Figure4_FAMILY_barplots_",group,".pdf",sep=""), onefile=FALSE)
dev.off()
}



#Figure 6 - RDA -----
figure6 = multi_panel_figure(width = 14,height = 28, unit = "inch",rows = 4,columns = 2)

figure6 %<>% fill_panel("figures/fungi/Figure6fs_RDA_Tomato.pdf",row = 1, column = 1)
figure6 %<>% fill_panel("figures/fungi/Figure6fr_RDA_Tomato.pdf",row = 2, column = 1)
figure6 %<>% fill_panel("figures/bacteria/Figure6bs_RDA_Tomato.pdf",row = 3, column = 1)
figure6 %<>% fill_panel("figures/bacteria/Figure6br_RDA_Tomato.pdf",row = 4, column = 1)

figure6 %<>% fill_panel("figures/fungi/Figure6fs_RDA_Pepper.pdf",row = 1, column = 2)
figure6 %<>% fill_panel("figures/fungi/Figure6fr_RDA_Pepper.pdf",row = 2, column = 2)
figure6 %<>% fill_panel("figures/bacteria/Figure6bs_RDA_Pepper.pdf",row = 3, column = 2)
figure6 %<>% fill_panel("figures/bacteria/Figure6br_RDA_Pepper.pdf",row = 4, column = 2)

figure6 %>% save_multi_panel_figure(filename = "figures/Figure6_rda.pdf")


#Figure5 - alpha -----
figure5 = multi_panel_figure(width = 18,height = 14, unit = "inch",rows = 2,columns = 2)

figure5 %<>% fill_panel("figures/fungi/Figure4fs_alpha.pdf",row = 1, column = 1)
figure5 %<>% fill_panel("figures/fungi/Figure4fr_alpha.pdf",row = 1, column = 2)
figure5 %<>% fill_panel("figures/bacteria/Figure4bs_alpha.pdf",row = 2, column = 1)
figure5 %<>% fill_panel("figures/bacteria/Figure4br_alpha.pdf",row = 2, column = 2)
figure5 %>% save_multi_panel_figure(filename = "figures/Figure5_alpha.pdf")


#Tree of candidate ASVs with taxonomy infos---
conditions = c("br","bs","fr","fs")
candidate.ASV.all = NULL
for(i in 1:4)
{
  candidate.ASV = read.table(paste("results/candidate.ASV.",conditions[i],".txt",sep = ""),stringsAsFactors = F)
  if(ncol(candidate.ASV) == 11) {candidate.ASV = candidate.ASV[,c(1:11,11)];candidate.ASV[,11] = "na";colnames(candidate.ASV)[11:12]  = c("Species","X")}
    candidate.ASV$condition = conditions[i]
  print(dim(candidate.ASV))
  candidate.ASV.all = rbind(candidate.ASV.all,candidate.ASV)
}

#collapse taxo info.
candidate.ASV.all$names_collapse = NULL
for(i in 1:80)
  {
    candidate.ASV.all$names_collapse[i]=paste(candidate.ASV.all[i,c(1,6:11)],collapse = ";")
    candidate.ASV.all$names_collapse[i]=gsub(";NA","",candidate.ASV.all$names_collapse[i],ignore.case=T)
    
    #how far along the taxo did you get...    
    len = length(gregexpr(";",candidate.ASV.all$names_collapse[i])[[1]])
    len.exact = paste("(",colnames(candidate.ASV.all)[5+len],")",sep= "")
    
    #get only the ASV + last part of the taxonomy
    temp_taxo=strsplit(candidate.ASV.all$names_collapse[i],split = ";")[[1]]
    if(length(temp_taxo)==1) candidate.ASV.all$names_collapse[i] = temp_taxo
    if(length(temp_taxo)>1) candidate.ASV.all$names_collapse[i] = paste(temp_taxo[1],temp_taxo[length(temp_taxo)],len.exact,collapse = ";")
    }

#plot trees
for(taxa in c("Bacteria","Fungi"))
{
  ###get phylo ----
  seqs <- candidate.ASV.all[candidate.ASV.all$Kingdom==taxa,12]
  names(seqs) <- candidate.ASV.all$names_collapse[candidate.ASV.all$Kingdom==taxa] # This propagates to the tip labels of the tree
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
  phang.align <- phyDat(as(alignment, "matrix"), type="DNA") #sequence alignments
  dm <- dist.ml(phang.align) #build distance matrix
  treeNJ <- NJ(dm) # Note, tip order != sequence order

  #plot phylo based on distance matrix----
  tip.color = c(1:80)
  for(i in 1:length(tip.color))
  {
    #get some nice colors to show up.
    if((candidate.ASV.all[i,4] == "Tomato") & (candidate.ASV.all[i,13] == "br")) tip.color[i] = "darkred"
    if((candidate.ASV.all[i,4] == "Tomato") & (candidate.ASV.all[i,13] == "bs")) tip.color[i] = "darkorange2"
    if((candidate.ASV.all[i,4] == "Pepper") & (candidate.ASV.all[i,13] == "br")) tip.color[i] = "cyan4"
    if((candidate.ASV.all[i,4] == "Pepper") & (candidate.ASV.all[i,13] == "bs")) tip.color[i] = "darkblue"
  
    if((candidate.ASV.all[i,4] == "Tomato") & (candidate.ASV.all[i,13] == "fr")) tip.color[i] = "darkred"
    if((candidate.ASV.all[i,4] == "Tomato") & (candidate.ASV.all[i,13] == "fs")) tip.color[i] = "darkorange2"
    if((candidate.ASV.all[i,4] == "Pepper") & (candidate.ASV.all[i,13] == "fr")) tip.color[i] = "cyan4"
    if((candidate.ASV.all[i,4] == "Pepper") & (candidate.ASV.all[i,13] == "fs")) tip.color[i] = "darkblue"
  }
  #plot.phylo
  dev.new()
  plot.phylo(treeNJ,cex = 0.5,tip.color = tip.color[candidate.ASV.all$Kingdom==taxa],font = 2,main = taxa,xpd = T)
  if(taxa == "Fungi") legend(0.02,40,legend =  c("Tomato - root","Tomato - soil","Pepper - root","Pepper - soil"), fill = c("darkred","darkorange2","cyan4","darkblue"),bg = "#FFFFFF99")
  if(taxa == "Bacteria") legend(0.02,40,legend =  c("Tomato - root","Tomato - soil","Pepper - root","Pepper - soil"), fill = c("darkred","darkorange2","cyan4","darkblue"),bg = "#FFFFFF99")
  dev.print(device=pdf, paste("figures/",tolower(taxa),"/Figure7_",tolower(taxa),"_tree.pdf",sep = ""), onefile=FALSE)
  dev.off()
  }


#Figure 7 - trees -----
figure7 = multi_panel_figure(width = 14,height = 7, unit = "inch",rows = 1,columns = 2)
figure7 %<>% fill_panel("figures/bacteria/Figure7_bacteria_tree.pdf",row = 1, column = 2)
figure7 %<>% fill_panel("figures/fungi/Figure7_fungi_tree.pdf",row = 1, column = 1)
figure7 %>% save_multi_panel_figure(filename = "figures/Figure7_candidateASVs.pdf")

#Figure 7 - candidates
#fungi root
#candidates = NULL
#candidates = rbind(candidates,taxo.abundants[colnames(asv.filt.abundants.norm) == "ASV70",])
#candidates = rbind(candidates,taxo.abundants[colnames(asv.filt.abundants.norm) == "ASV191",])
#fungi - soil
#candidates = rbind(candidates,taxo.abundants[colnames(asv.filt.abundants.norm) == "ASV243",])
#candidates = rbind(candidates,taxo.abundants[colnames(asv.filt.abundants.norm) == "ASV264",])
#candidates = rbind(candidates,taxo.abundants[colnames(asv.filt.abundants.norm) == "ASV252",])

#fungi - soil
taxo[c(132,153),]

#fungi - root
taxo[c(17,19),]
#





