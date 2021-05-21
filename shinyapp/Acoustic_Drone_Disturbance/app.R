library(shiny)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(latex2exp)
library(dplyr)

incoherentSum <- function(levels){
  10*log10(sum(10^(levels/10)))
}

#audio - data frame of SPL levels in 10hz bands, first column is frequencies, then magnitudes 
#bins - data frame of the bins, with columns for the lower, center, and upper values
#Returns data frame of the center frequencies and binned dB values 
binAudioFreqs <- function(audio, bins, type = 'SPL'){
  freqs <- audio[,1]
  mags <- audio[,2]
  
  binnedAudio <- lapply(1:nrow(bins),function(x){
    freqsInBin <- which(freqs>=bins[x,1]&freqs<bins[x,3])
    n <- length(freqsInBin)
    bandwidth <- bins[x,3]-bins[x,1]
    #For SPL
    if (type == 'SPL'){
      incoherentSum(mags[freqsInBin])-10*log10(n*10)+10*log10(bandwidth)
    }
    #For PSD
    else if (type == 'PSD'){
      incoherentSum(mags[freqsInBin])-10*log10(n*10)
    }
  })
  binnedAudio <- do.call(rbind,binnedAudio)
  out <- data.frame(bins[,2],binnedAudio)
  names(out) <- c('Frequency','Magnitude')
  out
}

drone_recordings <- read.csv('Data/Clean/drone_recordings_min.csv')
ambient <- read.csv('Data/Clean/ambient_recordings_min.csv')

noise_floor <- read.csv('Data/Reference/microphone_noise_floor.csv')
third_octaves <- read.csv('Data/Reference/third_octave.csv')
third_octaves_spotted_seal <- read.csv('Data/Reference/third_octave_spottedSeal.csv')  #Special case where octaves are not aligned around 1kHz

drone_list <- unique(drone_recordings$Drone)

#Only use third octaves with center freq geq 100
third_octaves <- third_octaves[third_octaves$Center >= 100,]

#thirdOctaves - data frame with third octave frequencies and magnitudes
third2octaves <- function(thirdOctaves){
  if (nrow(thirdOctaves) %% 3 != 0)
    stop('Number of third octave bands must be a multiple of three')
  octaveDF <- thirdOctaves[seq(2,nrow(thirdOctaves),3),]
  for (i in 1:nrow(octaveDF)){
    octaveDF[i,2:ncol(octaveDF)] <- incoherentSum(thirdOctaves[((i-1)*3+1):((i-1)*3+3),2])
  }
  names(octaveDF) <- c('Frequency','Magnitude')
  octaveDF
}

#soundCurve - data frame with frequencies and magnitudes, specify band type in 'format'  
#weights - data frame of frequencies (in Hz) and weights (dB)
#format - '10','Third', or 'Octave'
#bins - data frame of the bins, with columns for the lower, center, and upper values (only used if format=10) 
calcDBweighted <- function(soundCurve, weights, format, bins=third_octaves){
  if (nrow(weights) == 0 | ncol(weights) != 2)
    stop('Check that the weights are an nx2 data frame!')
  if (sum(weights[,2])>0)
    stop('Check that the weights are negative!')
  if (!format %in% c('Third','Octave','10')){
    stop('Only accepted binning is Third, Octave, or 10')
  }
  if (format != 'Octave'){
    if (format == '10'){
      soundCurve <- binAudioFreqs(soundCurve,bins)
    }
    soundCurve <- third2octaves(soundCurve)
  }
  weights <- weights[weights[,1] %in% soundCurve[,1],]
  soundCurve <- soundCurve[soundCurve[,1] %in% weights[,1],]
  if (nrow(weights)==0 | nrow(soundCurve)==0)
    stop('Check that the frequencies match between the sound curve and weights!')
  weightedCurve <- soundCurve[,2]+weights[,2]
  weightedCurve[weightedCurve < 0] <- 0 #Negative implies below hearing threshold so zero it out
  incoherentSum(weightedCurve)
}

#Calculate weighted dB vs distance in a plane for multiple drones
calcDBvsDistance <- function(drones, plane, weights, bins=third_octaves){
  all_drones <- lapply(drones, function(x){
    recordings <- filter(drone_recordings, Drone == x, Plane == plane)
    distances <- unique(recordings$Distance)
    one_drone <- lapply(distances, function(y){
      one_flight <- filter(recordings, Distance == y)
      weighted_spl <- calcDBweighted(one_flight,weights,'10',bins)
      out <- data.frame(x,y,weighted_spl)
      names(out) <- c('Drone','Distance','SPL')
      out
    })
    one_drone <- do.call(rbind, one_drone)
  })
  do.call(rbind, all_drones)
}

myPalette <- c('#000000',brewer.pal(12,'Paired')[c(2,4,6,7,8,10,12)])
plotDBvsDistance <- function(drones, species, bins=third_octaves){
  weights <- getSpeciesWeights(species)
  one_species <- calcDBvsDistance(drones, 'Vertical', weights, bins=bins)
  weighted_ambient <- calcDBweighted(ambient, weights,'10',bins)
  ambient_append <- data.frame(c('Ambient','Ambient'),
                               c(min(one_species$Distance),max(one_species$Distance)),
                               c(weighted_ambient,weighted_ambient))
  names(ambient_append) <- names(one_species)
  one_species <- rbind(one_species, ambient_append)
  one_species$Drone <- factor(one_species$Drone, levels = c('Ambient',drones)) #Legend order
  ggplot(data=one_species, aes(x=Distance,y=SPL,group=Drone)) +
    geom_point(aes(color=Drone, size=Drone)) +
    geom_line(aes(color=Drone,linetype=Drone)) +
    #ylim(c(20,70.5)) +
    scale_size_manual(values=c(-1,rep(1.5,length(drones)))) +
    scale_linetype_manual(values=c("dashed", rep("solid",length(drones))))+
    scale_color_manual(values = myPalette[c(1,(which(drone_list %in% drones)+1))]) +
    scale_x_continuous(trans='log2', breaks = c(0,5,10,20,40,80,120)) +
    xlab('Altitude [m]') +
    ylab(paste0('dB(',species,')')) +
    ggtitle(species) +
    theme_bw() +
    theme(legend.title=element_blank()) + 
    guides(color = guide_legend(override.aes = list(shape = c(NA,rep(16,length(drones)))))) +
    theme(text = element_text(size = 16)) +
    theme(plot.title = element_text(size = 20)) 
}

#Convert the species audiogram values into weights to add to the sound spectra
getSpeciesWeights <- function(species_name){
  extractSpecies <- audiograms[audiograms$Common.name==species_name,4:5]
  extractSpecies[,1] <- extractSpecies[,1]*1000
  extractSpecies[,2] <- extractSpecies[,2]*-1
  names(extractSpecies)[1] <- "Hz"
  extractSpecies
}

audiograms <- read.csv('Data/Clean/audiograms_inspected.csv')
species <- as.character(unique(audiograms$Common.name))




#drone - the name of the drone as a character
#weights - data frame of frequencies (in Hz) and weights (dB)
#bins - data frame of the bins, with columns for the lower, center, and upper values 
#threshold_slope - the slope (in -dB/m) that the best fit line must be less steep than
#            to qualify as acceptable altitude
#threshold_spl - dB threshold that is automatically acceptable
#loudness - set as TRUE if you want to output the dB(species) at the recommended altitude
calcFlyingAltitude <- function(drone, weights, bins=third_octaves, threshold_slope=0.01, threshold_spl=40, loudness=FALSE){
  sound_curve <- calcDBvsDistance(drone,'Vertical',weights,bins)
  
  flat_alt <- max(sound_curve$Distance) #Default advisable altitude at max altitude
  for (i in 1:(nrow(sound_curve)-1)){
    #Remove this condidtional to produce the table without the SPL threshold (Supporting Information) 
    if (sound_curve$SPL[i] <= threshold_spl){
      if (loudness == T)
        return(sound_curve$SPL[i])
      else
        return(sound_curve$Distance[i])
    }
    slope <- lm(sound_curve$SPL[i:nrow(sound_curve)]~sound_curve$Distance[i:nrow(sound_curve)])$'coefficients'[2]
    if (slope >= -threshold_slope){
      flat_alt <- sound_curve$Distance[i]
      break
    }
  }
  if (loudness == T)
    sound_curve$SPL[i]
  else
    flat_alt
}

makeAdvisableAltitudeTable <- function(drones, species, bins=third_octaves){
  altitudes <- data.frame(lapply(drones,calcFlyingAltitude, getSpeciesWeights(species), bins))
  names(altitudes) <- drones
  spls <- data.frame(lapply(drones,calcFlyingAltitude, getSpeciesWeights(species), bins, loudness=T))
  names(spls) <- drones
  full_table <- rbind(altitudes,spls)
  rownames(full_table) <- c('Advisable Altitude [m]','Weighted SPL [dB(species)]')
  full_table
}






# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Advisable altitudes for UAV flights over species"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      checkboxGroupInput("drones", 
                         'UAV Models',
                         choiceNames = drone_list,
                         choiceValues = drone_list),
      
      
      selectInput("species",
                  label = "Choose a species",
                  choices = species, selected='Reindeer')
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      tableOutput("advisable_altitude"),
      plotOutput(outputId = "species_plot")
      
    )
  )
)


server <- function(input, output) {
  output$advisable_altitude <- renderTable({ 
    if (length(input$drones) > 0) {
      if (input$species == 'Spotted Seal'){
        makeAdvisableAltitudeTable(input$drones, input$species, third_octaves_spotted_seal)
      } else{
        makeAdvisableAltitudeTable(input$drones, input$species)
      }
    }
  }, rownames=T, striped=T)
  output$species_plot <- renderPlot({
    if (length(input$drones) > 0) {
      if (input$species == 'Spotted Seal'){
        plotDBvsDistance(input$drones, input$species, third_octaves_spotted_seal)
      } else{
        plotDBvsDistance(input$drones, input$species)
      }
    }
  })
  
}

shinyApp(ui = ui, server = server)