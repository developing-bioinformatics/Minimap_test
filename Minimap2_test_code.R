library(taxonomizr)
library(thacklr)
library(dplyr)
library(forcats)
library(ggplot2)
# run minimap2
#system('fastq-dump SRR11043478')
#system('minimap2 -t32 /usr/share/data/ncbi/nt/nt.fa SRR11043478.fastq > approx-mapping.paf')
file="approx_mapping.paf"
file="test.sam"
X=read_paf(file, max_tags=20)
data=read.table(file, header=FALSE, sep= '\t')
#target=X[1:100, 6]
target=X[, 6]
accid = as.character(target)
taxaNodes<-read.nodes.sql("/usr/share/data/taxonomizr/nodes.dmp")
taxaNames<-read.names.sql("/usr/share/data/taxonomizr/names.dmp")
ids<-accessionToTaxa(target$target_name, '/usr/share/data/taxonomizr/accessionTaxa.sql')
taxlist=getTaxonomy(ids, taxaNodes, taxaNames)
print(taxlist)
Y=cbind(X,taxlist)
print(Y)
topmap=filter(Y, map_quality==0)
topmap=filter(topmap,dv<.1)
print(topmap)
cltop = topmap %>% 
  group_by(query_name) %>% 
  top_n(n=1, wt=map_length) %>%
  (ggplot(data=cltop) +
     geom_bar(aes(x=fct_infreq(genus))) +
     theme_minimal() +
     theme(    
       axis.text.x  = element_text(angle = 45, hjust=1)
     ) +
     xlab('')
   
  )