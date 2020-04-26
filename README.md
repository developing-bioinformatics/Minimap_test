#Gitclone this repository 
https://github.com/developing-bioinformatics/Minimap_test.git

#BLAST analysis

#First a BLAST database must be created, in this project we already had access to one. 

bl <- blast(db="/projectnb/ct-shbioinf/blast/nt.fa") 

#This line would need to be edited with their own nucleotide database for anyone who wishes to use it

#After that all that is to change this line inside rBLAST.R
srr=c('SRR11043497')
#to whatever eDNA sample file you want to sequence where it says 'SRR...'
#Then if your using a HPC cluster you can just use the BLAST_SubmitScript to run rBLAST.R otherwise you just run the entire file from there in your console.

#A singular png should be created with the top hits for families and genus.



#Minimap2 analysis 

# to create the index for Minimap2 also requires you to gitclone the Minimap2 repository before proceeding 
git clone https://github.com/lh3/minimap2

#Using Index_submit will create an index in the directory specified 
# make sure to alter the line "SRR11043497.fastq" with the appropriate sample file you want to sequence.

#after the index is completed as well as the paf file parsed by Minimap2, you can utilize the file Minimap.R to analyze and produce a png of the top hits for family and genus from your sample. 

