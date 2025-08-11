#This is for aligning any number of genes using ClustalW, ClustalOmega, or Muscle algorithms

#import aligner
library(msa)

#import file containing genes. Files need rows to represent ecotype and columns to represent genes
list <- read.csv(file = "C:/Users/kileegza/Documents/Non-Github analysis/Ecotype Variant Data/Alignments/Ecotype alignment_fasta.csv", row.names = 1, header = TRUE)

#simple for loop to go through the 231 columns present in the data set. Matrix is subset into Gene, which represents the gene
#to be aligned across the different ecotypes. 

  for (i in 1:231){
    
    Gene = as.character(list[1:9, i])
    
   
    
    out = msa(Gene, type ="dna")     #run the alignment using default settings with alignment type set to dna
    
    rownames(out) = paste(row.names(list), colnames(list[i]))    #set rownames equal to ecotype name followed by gene name
    
    #write the file to a folder. Name represents whichever gene was aligned across the ecotypes
    writeXStringSet(unmasked(out), file= paste("C:/Users/kileegza/Documents/Non-Github analysis/Ecotype Variant Data/Alignments/", colnames(list[i]), ".txt"))
    
  }
 
  



  



