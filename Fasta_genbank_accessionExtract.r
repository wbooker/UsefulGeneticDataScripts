Files <- c("Maturase", "rbcl", "ycf2", "atpb", "ndhf", "psba", "trnl", "ITS1")
mytree <- read.nexus("C:/Users/Will Booker/Documents/Macroevolution/AllFinal.nex")
list <- c("Species",mytree$tip.label)
GenesUsed <- c("Species","Maturase", "rbcl", "ycf2", "atpb", "ndhf", "psba", "trnl", "ITS1")
AllMat <- matrix(nrow=61, ncol=9)
AllMat[,1] <- list
AllMat[1,] <- GenesUsed
for (j in 1:length(Files)){
  gene<- Files[j]
  filePath <- paste(c("C:/Users/Will Booker/Documents/Macroevolution/Dataset/",gene,".fasta"), collapse = "")
  s = readDNAStringSet(filePath)
  names <- s@ranges@NAMES
  newmat <- matrix(nrow = length(names), ncol=2)
  for (i in 1:length(names)){
    temp1<- strsplit(names[i], " ")
    temp2<-temp1[[1]][1]
    temp3 <- strsplit(temp2, "|", fixed = TRUE)
    Accession <- temp3[[1]][4]
    SpeciesName <- paste(temp1[[1]][2], temp1[[1]][3], sep = "_")
    rowRef <- which(grepl(SpeciesName, AllMat[,1]))
    colRef <- which(grepl(gene, AllMat[1,]))
    AllMat[rowRef,colRef] <- Accession
  
  }


}

for (j in 1:length(Files)){
  gene<- Files[j]
  filePath <- paste(c("C:/Users/Will Booker/Documents/Macroevolution/Dataset/Outs/",gene,".fasta"), collapse = "")
  s = readDNAStringSet(filePath)
  names <- s@ranges@NAMES
  newmat <- matrix(nrow = length(names), ncol=2)
  for (i in 1:length(names)){
    temp1<- strsplit(names[i], " ")
    temp2<-temp1[[1]][1]
    temp3 <- strsplit(temp2, "|", fixed = TRUE)
    Accession <- temp3[[1]][4]
    SpeciesName <- paste(temp1[[1]][2], temp1[[1]][3], sep = "_")
    rowRef <- which(grepl(SpeciesName, AllMat[,1]))
    colRef <- which(grepl(gene, AllMat[1,]))
    AllMat[rowRef,colRef] <- Accession
    
  }
  
  
}

write.csv(AllMat,"C:/Users/Will Booker/Documents/Macroevolution/Accession.csv")
