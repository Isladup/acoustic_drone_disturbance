library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(latex2exp)
library(dplyr)

setwd('C:/Users/Marcus/Documents/MPhil/Acoustic_Drone_Disturbance/Data/Repository/')

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

################################################################################
#PSD PLOTS (VERTICAL, HORIZONTAL) FOR ALL DRONES AT FOUR DISTANCES

#Plot magnitude vs frequency for all drones at a given distance/direction
plotPSD <- function(plane, distance){
  #Filter the recordings
  drone_plot <- drone_recordings %>% 
    filter(Plane == plane, Distance == distance) %>% 
    select('Frequency','Magnitude','Drone')
  
  #PSD for drone, ambient, and noise floor
  drones <- unique(drone_plot$Drone)
  drone_psd <- lapply(drones, function(x){
    one_flight <- filter(drone_plot, Drone == x)
    one_flight_third <- binAudioFreqs(one_flight, third_octaves, 'PSD')
    one_flight_third$Drone <- x
    one_flight_third
  })
  drone_psd <- do.call(rbind,drone_psd)
  ambient_psd <- binAudioFreqs(ambient,third_octaves,'PSD')
  ambient_psd$Drone <- 'Ambient'
  noise_floor_psd <- binAudioFreqs(noise_floor_df,third_octaves,'PSD')
  noise_floor_psd$Drone <- 'Noise Floor'
  plot_psd <- rbind(drone_psd, ambient_psd, noise_floor_psd)
  
  #Palette for plotting
  myPalette <- brewer.pal(12,'Paired')
  myPalette <- myPalette[c(2,4,6,7,8,10,12)]
  myPalette <- c('#696969','#000000',myPalette)
  
  plot_psd$Drone <- factor(plot_psd$Drone, levels = c('Noise Floor','Ambient',drones)) #Legend order
  
  out <- ggplot(data=plot_psd, aes(x=Frequency,y=Magnitude,group=Drone))+
    geom_point(aes(color=Drone)) +
    geom_line(aes(color=Drone)) +
    scale_color_manual(values = myPalette) +
    ylim(c(-5,50)) +
    scale_x_continuous(trans='log2', breaks = c(250,1000,4000,16000)) +
    xlab('Frequency [Hz]') +
    ylab(TeX('\\overset{PSD}{$\\lbrack$dB re: 20$\\mu$ Pa$^2$/Hz$\\rbrack$}')) +
    ggtitle(paste0(distance, ' metres')) + 
    theme_bw() +
    theme(legend.title=element_blank()) 
}


drone_recordings <- read.csv('Clean/drone_recordings_min.csv')
ambient <- read.csv('Clean/ambient_recordings_min.csv')

noise_floor <- read.csv('Reference/microphone_noise_floor.csv')
third_octaves <- read.csv('Reference/third_octave.csv')
third_octaves_spotted_seal <- read.csv('Reference/third_octave_spottedSeal.csv')  #Special case where octaves are not aligned around 1kHz

#Only use third octaves with center freq geq 100
third_octaves <- third_octaves[third_octaves$Center >= 100,]

#Process the noise floor of the microphone, interpolating to match the recordings
interpX <- log10(ambient$Frequency)
noise_floor_interp <- approx(log10(noise_floor[,1]), noise_floor[,2], xout=interpX)
noise_floor_df <- data.frame(10^noise_floor_interp$x, noise_floor_interp$y)
names(noise_floor_df) <- c('Frequency','Magnitude')
noise_floor_df$Magnitude[2017:length(noise_floor_df$Magnitude)] <- noise_floor_df$Magnitude[2016]

h5 <- plotPSD('Horizontal',5)
h35 <- plotPSD('Horizontal',35)
h70 <- plotPSD('Horizontal',70)
h120 <- plotPSD('Horizontal',120)

v5 <- plotPSD('Vertical',5)
v35 <- plotPSD('Vertical',35)
v70 <- plotPSD('Vertical',70)
v120 <- plotPSD('Vertical',120)

figure <- ggarrange(v5,v35,v70,v120, nrow=2, ncol=2, common.legend=TRUE, legend='bottom')
annotate_figure(figure,
                top = text_grob("UAV Audio Recordings: Vertical", color = "black", face = "bold", size = 16)
)
figure <- ggarrange(h5,h35,h70,h120, nrow=2, ncol=2, common.legend=TRUE, legend='bottom')
annotate_figure(figure,
                top = text_grob("UAV Audio Recordings: Horizontal", color = "black", face = "bold", size = 16)
)


################################################################################
#CALCULATE AND PLOT SPECIES-WEIGHTED SPL

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

#Convert the species audiogram values into weights to add to the sound spectra
getSpeciesWeights <- function(species_name){
  extractSpecies <- audiograms[audiograms$Common.name==species_name,4:5]
  extractSpecies[,1] <- extractSpecies[,1]*1000
  extractSpecies[,2] <- extractSpecies[,2]*-1
  names(extractSpecies)[1] <- "Hz"
  extractSpecies
}

audiograms <- read.csv('Clean/audiograms_inspected.csv')
species <- as.character(unique(audiograms$Common.name))

#Select the species of interest
species_of_interest <- species[c(4,10,13)]
species_weights <- lapply(species_of_interest,getSpeciesWeights)

#A-weighted decibels
dBAfreq <- c(125,250,500,1000,2000,4000,8000,16000)
dBAweight <- c(-16.1,-8.6,-3.2,0,1.2,1,-1.1,-6.6)
dBA <- data.frame(dBAfreq, dBAweight)

#Append dB(A) weightings
species_weights <- c(list(dBA),species_weights)
species_names <- c('Human', species_of_interest)

#Palette for plotting
myPalette <- c('#000000',brewer.pal(12,'Paired')[c(2,4,6,7,8,10,12)])

#Plot sound curves for the species of interest in both planes
lapply(c('Horizontal','Vertical'), function(plane){
  one_direction <- lapply(1:length(species_weights), function(x){
    #Drone sound curves
    drones <- unique(drone_recordings$Drone)
    one_species <- calcDBvsDistance(drones,plane,species_weights[[x]],third_octaves)
    #Weighted ambient
    weighted_ambient <- calcDBweighted(ambient, species_weights[[x]],'10',third_octaves)
    ambient_append <- data.frame(c('Ambient','Ambient'),
                                 c(min(one_species$Distance),max(one_species$Distance)),
                                 c(weighted_ambient,weighted_ambient))
    names(ambient_append) <- names(one_species)
    one_species <- rbind(one_species, ambient_append)
    one_species$Drone <- factor(one_species$Drone, levels = c('Ambient',drones)) #Legend order
    direction_plot <- ggplot(data=one_species, aes(x=Distance,y=SPL,group=Drone)) +
      geom_point(aes(color=Drone, size=Drone)) +
      geom_line(aes(color=Drone,linetype=Drone)) +
      ylim(c(20,70.5)) +
      scale_size_manual(values=c(-1,rep(1.5,7))) +
      scale_linetype_manual(values=c("dashed", rep("solid",7)))+
      scale_color_manual(values = myPalette) +
      scale_x_continuous(trans='log2', breaks = c(0,5,10,20,40,80,120)) +
      ggtitle(paste0(species_names[x])) +
      theme_bw() +
      theme(legend.title=element_blank()) + 
      guides(color = guide_legend(override.aes = list(shape = c(NA,rep(16,7)))))
    if (species_names[x] == 'Human'){
      direction_plot <- direction_plot + ylab('dB(A)')
    } else{
      direction_plot <- direction_plot + ylab(paste0('dB(',species_names[x],')')) 
    }
    if (plane == 'Horizontal'){
      direction_plot <- direction_plot + xlab('Distance [m]')
    } else{
      direction_plot <- direction_plot + xlab('Altitude [m]')
    }
    direction_plot
  })
  figure <- ggarrange(one_direction[[1]],one_direction[[2]],
                      one_direction[[3]],one_direction[[4]],
                      nrow=2,ncol=2,common.legend=TRUE,legend='bottom')
  annotate_figure(figure, top = text_grob(paste0("Species-Filtered UAV Sound Level: ",plane),
                                          color = "black", face = "bold", size = 16)
  )
})

################################################################################
#CALCULATE ADVISABLE ALTITUDES BASED ON VERTICAL SPECIES-WEIGHTED SOUND CURVES

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

#Calculate the advisable altitude (or loudness at that altitude) for all drones for each species
#Note special treatment for spotted.seal
drones <- unique(drone_recordings$Drone)
all_species_altitudes <- lapply(species, function(x){
  species_weights <- getSpeciesWeights(x)
  species_all_drones <- lapply(drones, function(y){
    if (x == 'Spotted Seal'){
      #calcFlyingAltitude(y,species_weights,third_octaves_spotted_seal,0.01,40)
      sprintf('%0.1f',calcFlyingAltitude(y,species_weights,third_octaves_spotted_seal,0.01,40,T))
    } else{
      #calcFlyingAltitude(y,species_weights,third_octaves,0.01,40)
      sprintf('%0.1f',calcFlyingAltitude(y,species_weights,third_octaves,0.01,40,T))
    }
  })
  species_all_drones <- data.frame(do.call(cbind,species_all_drones))
  names(species_all_drones) <- drones
  species_all_drones$Species <- x
  species_all_drones
})
all_species_altitudes <- do.call(rbind,all_species_altitudes)
#Make species the first column
all_species_altitudes <- all_species_altitudes[,c(ncol(all_species_altitudes),1:(ncol(all_species_altitudes)-1))]
write.csv(all_species_altitudes, 'species_advisable_altitudes.csv')
#write.csv(all_species_altitudes, 'species_advisable_altitudes_SPL.csv')

################################################################################
#CALCULATE SENSITIVITES TO THRESHOLDS (SUPPORTING INFORMATION)

sensitivityPlots <- lapply(drones, function(x){
  sensitivityEachDrone <- lapply(species_of_interest, function(y){
    #thresholdRange <- seq(0.005,0.015,0.001) #For varying slope threshold
    thresholdRange <- 35:45 #For varying SPL threshold
    sensitivityVaryThreshold <- lapply(thresholdRange,function(thresh){
      #calcFlyingAltitude(x,getSpeciesWeights(y),third_octaves,thresh,40) #Slope
      calcFlyingAltitude(x,getSpeciesWeights(y),third_octaves,0.01,thresh) #SPL
    })
    sensitivities <- do.call(rbind,sensitivityVaryThreshold)
    sensitivities <- data.frame(thresholdRange, sensitivities)
    names(sensitivities) <- c('Threshold','Altitude')
    sensitivities$animal <- y
    sensitivities
  })
  sensitivityEachDrone <- do.call(rbind,sensitivityEachDrone)
  #sensitivityEachDrone$Threshold <- sensitivityEachDrone$Threshold*5 #Slope
  ggplot(data=sensitivityEachDrone, aes(x=Threshold,y=Altitude,group=animal)) +
    geom_point(aes(color=animal)) +
    geom_line(aes(color=animal)) +
    scale_color_brewer(palette="Set2") +
    #xlab('Slope Threshold [dB/5m]') + #Slope
    xlab('Noise Threshold [dB(species)]') + #SPL
    ylab('Advisable Altitude [m]') +
    ggtitle(x) +
    ylim(c(0,120)) +
    theme(legend.title=element_blank())
  
})

ggarrange(sensitivityPlots[[1]],sensitivityPlots[[2]],sensitivityPlots[[3]],sensitivityPlots[[5]],
          sensitivityPlots[[4]],sensitivityPlots[[6]],sensitivityPlots[[7]],
          nrow=2,ncol=4,common.legend=TRUE,legend='bottom')

