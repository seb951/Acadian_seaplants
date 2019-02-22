

setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants")

#combining figures
library(RColorBrewer)
library(vegan)
#library(multipanelfigure)
library(magrittr)
library(dplyr)
library(DECIPHER) #this is a bioconductor packages
library(phangorn)
library(phyloseq)

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
  if(taxa == "Fungi") legend(0.55,42,legend =  c("Tomato - root","Tomato - soil","Pepper - root","Pepper - soil"),xpd = T, fill = c("darkred","darkorange2","cyan4","darkblue"),bg = "#FFFFFFCC")
  if(taxa == "Bacteria") legend(0.25,30,legend =  c("Tomato - root","Tomato - soil","Pepper - root","Pepper - soil"), fill = c("darkred","darkorange2","cyan4","darkblue"),bg = "#FFFFFFCC")
  dev.print(device=pdf, paste("figures/",tolower(taxa),"/Figure7_",tolower(taxa),"_tree.pdf",sep = ""), onefile=FALSE)
  dev.off()
  }

###save the unknown FUNGI candidates to BLAST them.
candidate.ASV.f = candidate.ASV.all[41:80,]

candidate.ASV.fasta = rbind(
candidate.ASV.f[candidate.ASV.f[,1] == "ASV67",c(1,12)],
candidate.ASV.f[candidate.ASV.f[,1] == "ASV10",c(1,12)],
candidate.ASV.f[candidate.ASV.f[,1] == "ASV18",c(1,12)],
candidate.ASV.f[candidate.ASV.f[,1] == "ASV17",c(1,12)],
candidate.ASV.f[candidate.ASV.f[,1] == "ASV132",c(1,12)],
candidate.ASV.f[candidate.ASV.f[,1] == "ASV19",c(1,12)])

candidate.ASV.fasta[,1] = paste0(">",candidate.ASV.fasta[,1])

write.table(c(t(candidate.ASV.fasta)),"results/unknown_fungi_candidate.fasta",row.names = F, col.names = F, quote =F)

#Figure 7 - trees -----
#figure7 = multi_panel_figure(width = 14,height = 7, unit = "inch",rows = 1,columns = 2)
#figure7 %<>% fill_panel("figures/fungi/Figure7_fungi_tree.pdf",row = 1, column = 1)
#figure7 %<>% fill_panel("figures/bacteria/Figure7_bacteria_tree.pdf",row = 1, column = 2)
#figure7 %>% save_multi_panel_figure(filename = "figures/Figure7_candidateASVs.pdf")





