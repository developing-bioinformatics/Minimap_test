library(taxonomizr)
library(thacklr)
library(tibble)
library(stringr)
library(dplyr)
library(forcats)
library(ggplot2)
# run minimap2
#system('fastq-dump SRR11043478')
#system('minimap2 -t32 /usr/share/data/ncbi/nt/nt.fa SRR11043478.fastq > approx-mapping.paf')
file="approx_mapping.paf"
file="test2.paf"
X=read_paf(file, max_tags=20)
#data=read.table(file, header=FALSE, sep= '\t')
#target=X[1:100, 6]
target=X[, 6]
accid = as.character(target)
taxaNodes<-read.nodes.sql("/projectnb/ct-shbioinf/taxonomy/data/nodes.dmp")
taxaNames<-read.names.sql("/projectnb/ct-shbioinf/taxonomy/data/names.dmp")
ids<-accessionToTaxa(target$target_name, '/projectnb/ct-shbioinf/taxonomy/data/accessionTaxa.sql')
taxlist=getTaxonomy(ids, taxaNodes, taxaNames)
print(taxlist)
Y=cbind(X,taxlist)
print(Y)
topmap=filter(Y, map_quality==0)
topmap=filter(topmap,dv<.1)
print(topmap)
cltop = topmap %>% 
  group_by(query_name) %>% 
  top_n(n=1, wt=map_length)%>%
  group_by(genus)%>%
  summarise(count=n())%>%
  filter(count>100)



  (ggplot(data=cltop) +
     geom_col(aes(x=genus, y=count)) +
      scale_y_log10()+
     theme_minimal() +
     theme(    
       axis.text.x  = element_text(angle = 45, hjust=1)
     ) +
     xlab('')
   
  )
source("LCA_Minimap2.R")

temp=lca(topmap, parallel=TRUE, nclus=4)

run = temp %>% 
  group_by(query_name) %>% 
  top_n(n=1, wt=map_length)%>%
  group_by(genus)%>%
  summarise(count=n())%>%
  filter(count>5)

(ggplot(data=run) +
    geom_col(aes(x=genus, y=count)) +
    scale_y_log10()+
    theme_minimal() +
    theme(    
      axis.text.x  = element_text(angle = 45, hjust=1)
    ) +
    xlab('')
  
)

familyRun=temp %>% 
  group_by(query_name) %>% 
  top_n(n=1, wt=map_length)%>%
  group_by(family)%>%
  summarise(count=n())%>%
  filter(count>5)


(ggplot(data=familyRun) +
    geom_col(aes(x=family, y=count)) +
    scale_y_log10()+
    theme_minimal() +
    theme(    
      axis.text.x  = element_text(angle = 45, hjust=1)
    ) +
    xlab('')
  
)