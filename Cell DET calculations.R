setwd("C:/Users/avining/Downloads")
tracks <- read.csv("route_data_subset_gridded.csv")
library(dplyr)
library(spatstat)


get_cell_DET <- function(x, minl){
  
  x = as.numeric(x)
  #Depending on the dataset it may be desirable to filter out consecutive visits 
  #to the same flower. See function below and delete '#' in the line below to use
  #x = filterout(Ldata = x) #already done
  
  #-----set up matrix resembling a recurrence plot, where a 1 indicates a repeat 
  #-----visit and 0 indicates the absence of a repeat.
  #if (length(unique(x)) < minl) return(NA)
  #if (length(x) <= 3*minl) return(NA)
  #if (length(unique(x)) == length(x)) return(0)
  det1 = matrix(cbind(rep(x,length(x))),nrow=length(x))
  tdet = t(det1)
  det = ((det1 - tdet) == 0) * 1
  
  #set the main diagonal equal to zero so it won't be included in the calculation
  
  diag(det) = 0
  
  #Use spatstat package to create a 'countour map' of the matrix,
  #which assigns all sets of contiguous 1's a unique number
  yi <- as.im(det)
  ycOut <- connected(yi, background = 0)
  yc <- ycOut$v
  
  #Depending on the dataset it may be desirable to filter out diagonals perpendicular to #the main diagonal. Code is provided for the 'removeperpdiag' function below.
  #Delete "#" from the line below to filter out perpendicular diagonals
  
  #yc = removeperpdiag(yc,minl)
  
  #Note: this code may take several minutes to run for very long sequences
  
  #---- filter out short repeats: a 'trapline' should include more unique resources
  #---- than the minimum cutoff (minl)
  
  #make an alternative DET matrix that contains the resource IDs
  det2 = matrix(rep(x,nrow(det)),nrow=nrow(det),byrow=TRUE)*det
  recCounts <- as.data.frame(table(det2)/2) #get number of times each cell occurs in the top corner of the matrix
  #make a dataframe with the number of times each resource appears in a diagonal
  listofseq = data.frame(group = yc[1:length(yc)], seq=det2[1:length(det2)])
  #how many unique resources are in the diagonal
  uniquevisits = rowSums((table(listofseq)>0)*1)
  #only count diagonals with at least 'minl' number of unique resources
  longenough = (uniquevisits >= minl)*table(yc)
  listof_longenoug_seq <- listofseq[listofseq$group %in% names(longenough[! longenough == 0]),]
  repCounts <- as.data.frame(table(listof_longenoug_seq$seq)/2)
  cellDET <- as.data.frame(table(x))
  colnames(recCounts) <- c("cell", "recursions")
  colnames(repCounts) <- c("cell", "repeats")
  colnames(cellDET) <- c("cell", "visits")
  cellDET <- merge(cellDET, recCounts, by = "cell", all = TRUE)
  cellDET <- merge(cellDET, repCounts, by = "cell", all = TRUE)
  cellDET$repeats[is.na(cellDET$repeats)] <- 0
  cellDET <- mutate(cellDET, DET = repeats/recursions)
  return(cellDET)
}

cellDETs <- vector("list", length = length(levels(tracks$ID)))
names(cellDETs) <- levels(tracks$ID)
for(i in 1:length(levels(tracks$ID))) {
  sequence <- rle(filter(tracks, ID == levels(tracks$ID)[i])$cell)$values
  cellDETs[[i]] <- get_cell_DET(sequence, minl = 3)
}




