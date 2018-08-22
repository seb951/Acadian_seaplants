

setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants")

#combining figures
library(multipanelfigure)
library(magrittr)
library(dplyr)
library(DECIPHER)
library(phangorn)
library(phyloseq)

#Figure 4 - barplots -----
figure4 = multi_panel_figure(width = 28,height = 8, unit = "inch",rows = 1,columns = 4)

figure4 %<>% fill_panel("figures/bacteria/Figure4br_FAMILY_ASVabundance.pdf",row = 1, column = 1)
figure4 %<>% fill_panel("figures/bacteria/Figure4bs_FAMILY_ASVabundance.pdf",row = 1, column = 2)
figure4 %<>% fill_panel("figures/fungi/Figure4fr_FAMILY_ASVabundance.pdf",row = 1, column = 3)
figure4 %<>% fill_panel("figures/fungi/Figure4fs_FAMILY_ASVabundance.pdf",row = 1, column = 4)

figure4 %>% save_multi_panel_figure(filename = "figures/Figure4_barplots.pdf")

#Figure 6 - RDA -----
figure6 = multi_panel_figure(width = 28,height = 14, unit = "inch",rows = 2,columns = 4)

figure6 %<>% fill_panel("figures/bacteria/Figure6bs_RDA_Tomato.pdf",row = 1, column = 1)
figure6 %<>% fill_panel("figures/bacteria/Figure6bs_RDA_Pepper.pdf",row = 2, column = 1)
figure6 %<>% fill_panel("figures/bacteria/Figure6br_RDA_Tomato.pdf",row = 1, column = 2)
figure6 %<>% fill_panel("figures/bacteria/Figure6br_RDA_Pepper.pdf",row = 2, column = 2)
figure6 %<>% fill_panel("figures/fungi/Figure6fs_RDA_Tomato.pdf",row = 1, column = 3)
figure6 %<>% fill_panel("figures/fungi/Figure6fs_RDA_Pepper.pdf",row = 2, column = 3)
figure6 %<>% fill_panel("figures/fungi/Figure6fr_RDA_Tomato.pdf",row = 1, column = 4)
figure6 %<>% fill_panel("figures/fungi/Figure6fr_RDA_Pepper.pdf",row = 2, column = 4)

figure6 %>% save_multi_panel_figure(filename = "figures/Figure6_rda.pdf")


#Figure5 - alpha -----
figure5 = multi_panel_figure(width = 18,height = 14, unit = "inch",rows = 2,columns = 2)
figure5 %<>% fill_panel("figures/bacteria/Figure4bs_alpha.pdf",row = 1, column = 1)
figure5 %<>% fill_panel("figures/bacteria/Figure4br_alpha.pdf",row = 1, column = 2)
figure5 %<>% fill_panel("figures/fungi/Figure4fs_alpha.pdf",row = 2, column = 1)
figure5 %<>% fill_panel("figures/fungi/Figure4fr_alpha.pdf",row = 2, column = 2)
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
  legend(max(treeNJ$edge.length),40,legend =  c("Tomato - root","Tomato - soil","Pepper - root","Pepper - soil"), fill = c("darkred","darkorange2","cyan4","darkblue"))
  dev.print(device=pdf, paste("figures/",tolower(taxa),"/Figure7_",tolower(taxa),"_tree.pdf",sep = ""), onefile=FALSE)
  dev.off()
  }
