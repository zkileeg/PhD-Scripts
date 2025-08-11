#Take list of genes in a fasta file and find which are the most similar to the gene/protein of interest

#load required
library(RecordLinkage)
library(stringdist)
library(stringr)
library(ggvis)




#Get sequence data from csv files. Change depending on requirements

protein <- read.csv("C:/Users/kileegza/Documents/Non-Github analysis/Ecotype Variant Data/lrrrlklist.csv", sep = ",", header = TRUE)

stringvalues <- read.csv("C:/Users/kileegza/Documents/Non-Github analysis/Ecotype Variant Data/Sha/Sha_protlist.csv", header = FALSE, sep = ",")

#Test file
#test <- read.csv("C:/Users/kileegza/Documents/Non-Github analysis/test.csv", header = TRUE, sep = ",")

#flatten list and coerce values into a character vector
compare = as.character(unlist(stringvalues$V1))


#function closest match. Take string list containing LRRs and, sequentially, calculate the 'string' distance 
#between the protein sequence and the 22000 sequences from the ecotype proteome. Return the string containing
#the highest value (most similar sequence)

Closestmatch = function(stringlist, string) {
  
  
  distance = levenshteinSim(stringlist, string)
  string[distance == max(distance)]
}



#sequentially go through protein codes and send it to Closest match to get string from ecotype closest matching columbia sequence
Output = lapply(as.character(protein$Protein), Closestmatch, string = compare)


#write output to a csv
write.csv(Output,"C:/Users/kileegza/Documents/Non-Github analysis/Ecotype Variant Data/Sha_prot_out.csv")
