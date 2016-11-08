root <- "C:/Users/Will Booker/Dropbox/Documents/Research/PhD/Migrate/"
lociTable <- as.matrix(read.csv(paste((c(root, "LociTable.csv")), collapse = ""), header = FALSE))
individualsTable <- as.matrix(read.csv(paste((c(root, "IndTable.csv")), collapse = ""), header = FALSE))
proj <- "T177"
npops <- 3

############Split by Locus
allLociMat <- matrix(nrow=1, ncol = 3)
LociLength <- matrix(nrow = 1, ncol = length(lociTable[,1]))
for (i in 1:length(lociTable[,1])){
  filepath <- paste((c(root,proj,"_L",as.character(lociTable[i,1]),".fasta")), collapse = "")
  #filepath <- "C:/Users/Will Booker/Dropbox/Documents/Research/PhD/Migrate/T177_L15.fasta"
  k <-scan(filepath, what = "complex")
  IndividualsMat <- matrix(ncol = 1, nrow=length <- (as.numeric(length(k))/2))
  allelesMat <- matrix(ncol = 1, nrow=length <- (as.numeric(length(k))/2))
  counter <- 1
  locusMat <- matrix(nrow=1, ncol=3)
  ############Split all alleles in Locus
  for (j in seq(1,(length(k)-1),2)){
    IndividualsMat[counter,1] <- substring(strsplit(k[j],"_", fixed = TRUE)[[1]][1],2)
    allelesMat[counter,1] <- k[j+1]
    counter <- counter+1
  }
    ###########Look at only individuals of Interest
    for(g in 1:length(individualsTable[,1])){
      alleleCount <- 1
      individualRef <- which(grepl(as.character(individualsTable[g,1]), IndividualsMat))
      #########Correct number of alleles
      for(z in 1:length(individualRef)){
        TempMat <- matrix(nrow = 1, ncol = 3)
        TempMat[1,1]<- paste(c(as.character(individualsTable[g,1]),"_",alleleCount), collapse = "")
        TempMat[1,2]<- allelesMat[as.numeric(individualRef[alleleCount]),1]
        TempMat[1,3]<- individualsTable[g,2]
        locusMat <- rbind(locusMat, TempMat)  
        alleleCount <- alleleCount + 1
      }
      
      alleleCount <- 1 
      
    }

  locusMat <- locusMat[2:length(locusMat[,1]),]
  LociLength[1,i] <- nchar(locusMat[1,2])
  allLociMat <- rbind(allLociMat,locusMat)
  ############Individual locu and individuals mat good, 
  ############just need to get it to go through all 
  ############loci and make a big file

}
lociStr <- paste(LociLength, sep = " ", collapse = " ")
allLociMat <- allLociMat[2:length(allLociMat[,1]),]
allLociMat[,1]<-sprintf(fmt = "%-10s",allLociMat[,1])

########## Need to change pops info here########
finalMat <- rbind(c(npops, length(lociTable[,1])), c(lociStr, ""),  c((length(allLociMat[allLociMat[,3]==1,1])/length(lociTable[,1])), "NE"), allLociMat[allLociMat[,3]==1,1:2], c((length(allLociMat[allLociMat[,3]==2,1])/length(lociTable[,1])), "NW"), allLociMat[allLociMat[,3]==2,1:2], c((length(allLociMat[allLociMat[,3]==3,1])/length(lociTable[,1])), "SW"), allLociMat[allLociMat[,3]==3,1:2])
##############################################

write.table(finalMat, file = paste(c(root,"test.txt"),collapse = ""), sep = "\t", row.names=FALSE, col.names=FALSE, quote = FALSE)

