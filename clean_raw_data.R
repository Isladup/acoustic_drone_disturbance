library(gtools)
library(dplyr)
library(stringr)

setwd('C:/Users/Marcus/Documents/MPhil/Acoustic_Drone_Disturbance/Data/Repository/')

#csv - The full filename of the CSV file containing the audio recording
#Returns data frame with columns for frequency and magnitude
readAudioCSV <- function(csv){
  findStart <- read.csv(csv)
  startInd <- which(rownames(findStart)=='Frequency')
  ambientWave <- read.csv(csv,skip=startInd)
  #Add 20.9 to simulate switch from low to medium gain
  ambientWave[,2] <- ambientWave[,2]+20.9
  ambientWave
}

################################################################################
#AMBIENT RECORDINGS

#Process the ambient recordings into a single, minimum ambient wave
ambients <- list.files('Raw/Wytham Drone Recordings/',pattern='mbient',recursive=TRUE, full.names=TRUE)

ambientList <- lapply(ambients, function(x){
  ambientMag <- readAudioCSV(x)$Magnitude
}) 
ambientAll <- do.call(cbind,ambientList)
ambientMin <- apply(ambientAll,1,min) #Min at each frequency 
ambientClean <- readAudioCSV(ambients[1])#Attach the frequency values
ambientClean$Frequency <- ambientClean$Frequency
ambientClean$Magnitude <- ambientMin
write.csv(ambientClean, 'Clean/ambient_recordings_min.csv', row.names = F)

################################################################################
#DRONE RECORDINGS

#Process each set of drone recordings into a single minimum spectrum

#List the recordings and remove ambients and the 0 distance
droneFolders <- list.files('Raw/Wytham Drone Recordings/') #separate folder for each drone
droneFiles <- mixedsort(list.files('Raw/Wytham Drone Recordings/',recursive=TRUE, full.names=TRUE))
droneFiles <- droneFiles[!grepl('mbient',droneFiles)]
droneFiles <- droneFiles[!grepl('/0',droneFiles)]

#Also remove Phantom at 100 and Mavic Mini at 40 metres because they're problematic
droneFiles <- droneFiles[!grepl('Phantom 4/Vertical/100', droneFiles)]
droneFiles <- droneFiles[!grepl('Mavic Mini/Vertical/40', droneFiles)]

#Read in the drone recordings and categorise by drone, plane, and distance
droneRecs <- lapply(droneFiles, function(x){
  recording <- readAudioCSV(x)
  #Drone
  recording$Drone <- unlist(lapply(droneFolders,function(y){
    if (grepl(paste0(y,'/'),x)==TRUE) #Include the / to avoid catching names within names
      y
  }))
  #Plane
  recording$Plane <- unlist(lapply(c('Horizontal','Vertical'),function(y){
    if (grepl(y,x)==TRUE)
      y
  }))
  #Distance
  startDist <- unlist(gregexpr(pattern ='al/',x))+3 #After horizontal/vertical
  i <- startDist
  while (TRUE) {
    check <- substr(x,i,i)
    if (check == '_' || check == '.' || check == 'm' || check == ' ')
      break #Break when reached not a number
    i <- i+1
  }
  recording$Distance <- substr(x,startDist,(i-1))
  recording
})
droneRecsAll <- do.call(rbind,droneRecs)

#For each flight, generate one spectrum with the minimum magnitude at each frequency
replicate_string <- paste(droneRecsAll$Drone,droneRecsAll$Plane,droneRecsAll$Distance,sep='_')
unique_flights <- unique(replicate_string)
min_recordings <- lapply(unique_flights, function(x){
  one_flight <- droneRecsAll[replicate_string==x,]
  test <- one_flight %>% 
    group_by(Frequency) %>% 
    summarise_all(min)
})
clean_recordings <- do.call(rbind,min_recordings)
write.csv(clean_recordings, 'Clean/drone_recordings_min.csv', row.names = F)

################################################################################
#AUDIOGRAMS

#Capitalise the first letter of every word
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

#Filter out the audiograms based on the binary flag done by inspection
audiograms <- read.csv('Raw/Audiograms/Audiograms_2021-01-27.csv')
audiograms <- audiograms[audiograms$X==1,]
#Format the common name
audiograms$Common.name <- unlist(lapply(str_replace_all(audiograms$Common.name,'\\.',' '),CapStr))

write.csv(audiograms, 'Clean/audiograms_inspected.csv', row.names = F)
