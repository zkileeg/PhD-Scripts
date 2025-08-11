#GC3 content calculator. GC3 = the guanine-cytosine content at the third position of a codon. Synonymous with the wobble 
#position. 
#The more recombination ,repair, and positive selection that occurs, the higher the GC3 content. 
#This could be useful to predict evolutionary events


#Get cDNA file from text

#function for counting GC content

#get packages for counting string occurrence
library(stringr)



#get data from CSV containing LRR cDNA data
genes <- read.csv("C:/Users/kileegza/Documents/Non-Github analysis/LRRs_GC3.csv", header = TRUE, sep = ",", row.names = 1) 



#take the genes data frame from before and apply the GC_Count function to each entry, returning a gene count into list numbers
#coerce into a character vector first, then apply function GC_Count
numbers <- lapply(as.character(genes$cDNA), GC_Count)


#Find the number of Gs or Cs at the 3rd position of each codon
GC_Count <- function(gene, increment = 3){
  
  #get value of the third character and put it into i. Start at the zeroth element, and continute by 3 until the length of gene
  i <- seq(0, nchar(gene), increment)
  
  #create subset by combining value of letter from i (start at i finish at i), and continute until end
  subset <- paste0(mapply(substr, gene, i, i), collapse="")
 
  #count the number of Cs or Gs present in the string and put it into a percentage of the total
  GC3 <- ((str_count(subset, "c") + str_count(subset, "g")) / (nchar(gene)/3)) * 100
  
  
  
  return(GC3)
  
  
}


#Combine the elements from genes and add the column from numbers with GC3 content into a new column. Coerce list into vector for writing to csv
final <- genes
final$GC3 <- as.vector(unlist(numbers))

write.csv(final, "C:/Users/kileegza/Documents/Non-Github analysis/GC3.csv")
