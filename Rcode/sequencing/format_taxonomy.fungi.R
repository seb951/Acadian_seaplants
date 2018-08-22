setwd("/Users/jerry/Documents/CSBQ/hijri/Acadian_seaplants/reference_material")


taxo.fungi = read.table("vtx_taxonomy.txt",stringsAsFactors = F, sep = "\t", header = T)

#taxo reformat
format.taxo = c(1:nrow(taxo.fungi))
for(i in 1:length(format.taxo))
{
  format.taxo[i]=paste(taxo.fungi[i,c(4:8)],collapse = ";")
  format.taxo[i]=paste(">Fungi",";",format.taxo[i],"_",taxo.fungi[i,2],";",sep = "")
}


#save fasta
write.table(c(rbind(format.taxo,taxo.fungi[,3])),"vtx_taxonomy.reformatted.fasta",row.names =F, col.names =F, quote = F)


###format Unite DB ### 

unite.fasta = read.table("unite_singleton/sh_general_release_dynamic_s_01.12.2017.fasta",stringsAsFactors = F, sep = "\t", header = F)

#taxo reformat  (>Fungi;Paraglomeromycetes;Paraglomerales;Paraglomeraceae;Paraglomus;sp._VTX00001;)
format.taxo = c(1:nrow(unite.fasta))
for(i in seq(1,length(format.taxo),by = 2))
{
  format.taxo[i]=paste(gsub(".__","",(strsplit(unite.fasta[i,],";")[[1]][2:7])),collapse = ";")
  
  format.taxo[i]=paste(">Fungi",";",format.taxo[i],";",sep = "")
}

unite.fasta[seq(1,length(format.taxo),by = 2),1] = format.taxo[seq(1,length(format.taxo),by = 2)]

#save fasta
write.table(unite.fasta,"unite_singleton/sh_general_release_dynamic_s_01.12.2017_reformated.fasta",row.names =F, col.names =F, quote = F)

