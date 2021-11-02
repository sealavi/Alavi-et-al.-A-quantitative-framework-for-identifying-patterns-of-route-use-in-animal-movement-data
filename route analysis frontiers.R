####Load data from Movebank

library(getPass)  
library(move)

pass <- getPass::getPass() ##keep password confidential
loginStored <- movebankLogin(username="Shauhin", password=pass)

Abby1 <- getMovebankData(study=1120749252 , animalName="Abby 4652", login=loginStored)
Abby1@data$individual_taxon_canonical_name=Abby1@idData$taxon_canonical_name
Abby1@data$individual_local_identifier=Abby1@idData$local_identifier
Abby1 <- data.frame(Abby1@data)

Abby2 <- getMovebankData(study=1120749252 , animalName="Abby 5767", login=loginStored)
Abby2@data$individual_taxon_canonical_name=Abby2@idData$taxon_canonical_name
Abby2@data$individual_local_identifier=Abby2@idData$local_identifier
Abby2 <- data.frame(Abby2@data)

Avery <- getMovebankData(study=1120749252 , animalName="Avery 4671", login=loginStored)
Avery@data$individual_taxon_canonical_name=Avery@idData$taxon_canonical_name
Avery@data$individual_local_identifier=Avery@idData$local_identifier
Avery <- data.frame(Avery@data)

BenBob <- getMovebankData(study=1120749252 , animalName="Ben Bob 4653", login=loginStored)
BenBob@data$individual_taxon_canonical_name=BenBob@idData$taxon_canonical_name
BenBob@data$individual_local_identifier=BenBob@idData$local_identifier
BenBob <- data.frame(BenBob@data)

Bob <- getMovebankData(study=1120749252 , animalName="Bob 4661", login=loginStored)
Bob@data$individual_taxon_canonical_name=Bob@idData$taxon_canonical_name
Bob@data$individual_local_identifier=Bob@idData$local_identifier
Bob <- data.frame(Bob@data)

Bonnie <- getMovebankData(study=1120749252 , animalName="Bonnie 4658", login=loginStored)
Bonnie@data$individual_taxon_canonical_name=Bonnie@idData$taxon_canonical_name
Bonnie@data$individual_local_identifier=Bonnie@idData$local_identifier
Bonnie <- data.frame(Bonnie@data)

Carlsberg <- getMovebankData(study=1120749252 , animalName="Carlsberg 4673", login=loginStored)
Carlsberg@data$individual_taxon_canonical_name=Carlsberg@idData$taxon_canonical_name
Carlsberg@data$individual_local_identifier=Carlsberg@idData$local_identifier
Carlsberg <- data.frame(Carlsberg@data)

Chibi <- getMovebankData(study=1120749252 , animalName="Chibi 4693", login=loginStored)
Chibi@data$individual_taxon_canonical_name=Chibi@idData$taxon_canonical_name
Chibi@data$individual_local_identifier=Chibi@idData$local_identifier
Chibi <- data.frame(Chibi@data)

Chloe <- getMovebankData(study=1120749252 , animalName="Chloe 4052", login=loginStored)
Chloe@data$individual_taxon_canonical_name=Chloe@idData$taxon_canonical_name
Chloe@data$individual_local_identifier=Chloe@idData$local_identifier
Chloe <- data.frame(Chloe@data)

Clementina <- getMovebankData(study=1120749252 , animalName="Clementina 4672", login=loginStored)
Clementina@data$individual_taxon_canonical_name=Clementina@idData$taxon_canonical_name
Clementina@data$individual_local_identifier=Clementina@idData$local_identifier
Clementina <- data.frame(Clementina@data)

Da_Vinci <- getMovebankData(study=1120749252 , animalName="Da Vinci 5764", login=loginStored)
Da_Vinci@data$individual_taxon_canonical_name=Da_Vinci@idData$taxon_canonical_name
Da_Vinci@data$individual_local_identifier=Da_Vinci@idData$local_identifier
Da_Vinci <- data.frame(Da_Vinci@data)

Eli <- getMovebankData(study=1120749252 , animalName="Eli 5765", login=loginStored)
Eli@data$individual_taxon_canonical_name=Eli@idData$taxon_canonical_name
Eli@data$individual_local_identifier=Eli@idData$local_identifier
Eli <- data.frame(Eli@data)

Ellie <- getMovebankData(study=1120749252 , animalName="Ellie 4668", login=loginStored)
Ellie@data$individual_taxon_canonical_name=Ellie@idData$taxon_canonical_name
Ellie@data$individual_local_identifier=Ellie@idData$local_identifier
Ellie <- data.frame(Ellie@data)

Emma <- getMovebankData(study=1120749252 , animalName="Emma 5762", login=loginStored)
Emma@data$individual_taxon_canonical_name=Emma@idData$taxon_canonical_name
Emma@data$individual_local_identifier=Emma@idData$local_identifier
Emma <- data.frame(Emma@data)

Fonta_Flora <- getMovebankData(study=1120749252 , animalName="Fonta Flora 4689", login=loginStored)
Fonta_Flora@data$individual_taxon_canonical_name=Fonta_Flora@idData$taxon_canonical_name
Fonta_Flora@data$individual_local_identifier=Fonta_Flora@idData$local_identifier
Fonta_Flora <- data.frame(Fonta_Flora@data)

Galena <- getMovebankData(study=1120749252 , animalName="Galena 5775", login=loginStored)
Galena@data$individual_taxon_canonical_name=Galena@idData$taxon_canonical_name
Galena@data$individual_local_identifier=Galena@idData$local_identifier
Galena <- data.frame(Galena@data)

Gamer <- getMovebankData(study=1120749252 , animalName="Gamer 5772", login=loginStored)
Gamer@data$individual_taxon_canonical_name=Gamer@idData$taxon_canonical_name
Gamer@data$individual_local_identifier=Gamer@idData$local_identifier
Gamer <- data.frame(Gamer@data)

Gillian <- getMovebankData(study=1120749252 , animalName="Gillian 4671", login=loginStored)
Gillian@data$individual_taxon_canonical_name=Gillian@idData$taxon_canonical_name
Gillian@data$individual_local_identifier=Gillian@idData$local_identifier
Gillian <- data.frame(Gillian@data)

Greg <- getMovebankData(study=1120749252 , animalName="Greg 4689", login=loginStored)
Greg@data$individual_taxon_canonical_name=Greg@idData$taxon_canonical_name
Greg@data$individual_local_identifier=Greg@idData$local_identifier
Greg <- data.frame(Greg@data)

Golliath <- getMovebankData(study=1120749252 , animalName="Golliath 5214", login=loginStored)
Golliath@data$individual_taxon_canonical_name=Golliath@idData$taxon_canonical_name
Golliath@data$individual_local_identifier=Golliath@idData$local_identifier
Golliath <- data.frame(Golliath@data)

Ibeth <- getMovebankData(study=1120749252 , animalName="Ibeth 4654", login=loginStored)
Ibeth@data$individual_taxon_canonical_name=Ibeth@idData$taxon_canonical_name
Ibeth@data$individual_local_identifier=Ibeth@idData$local_identifier
Ibeth <- data.frame(Ibeth@data)

Inez <- getMovebankData(study=1120749252 , animalName="Inez 5213", login=loginStored)
Inez@data$individual_taxon_canonical_name=Inez@idData$taxon_canonical_name
Inez@data$individual_local_identifier=Inez@idData$local_identifier
Inez <- data.frame(Inez@data)

Jeff <- getMovebankData(study=1120749252 , animalName="Jeff 5769", login=loginStored)
Jeff@data$individual_taxon_canonical_name=Jeff@idData$taxon_canonical_name
Jeff@data$individual_local_identifier=Jeff@idData$local_identifier
Jeff <- data.frame(Jeff@data)

Judy <- getMovebankData(study=1120749252 , animalName="Judy 4656", login=loginStored)
Judy@data$individual_taxon_canonical_name=Judy@idData$taxon_canonical_name
Judy@data$individual_local_identifier=Judy@idData$local_identifier
Judy <- data.frame(Judy@data)

Kyle <- getMovebankData(study=1120749252 , animalName="Kyle 4692", login=loginStored)
Kyle@data$individual_taxon_canonical_name=Kyle@idData$taxon_canonical_name
Kyle@data$individual_local_identifier=Kyle@idData$local_identifier
Kyle <- data.frame(Kyle@data)

Limon <- getMovebankData(study=1120749252 , animalName="Limon 5215", login=loginStored)
Limon@data$individual_taxon_canonical_name=Limon@idData$taxon_canonical_name
Limon@data$individual_local_identifier=Limon@idData$local_identifier
Limon <- data.frame(Limon@data)

Mario <- getMovebankData(study=1120749252 , animalName="Mario 5768", login=loginStored)
Mario@data$individual_taxon_canonical_name=Mario@idData$taxon_canonical_name
Mario@data$individual_local_identifier=Mario@idData$local_identifier
Mario <- data.frame(Mario@data)

Martinelli <- getMovebankData(study=1120749252 , animalName="Martinelli 5763", login=loginStored)
Martinelli@data$individual_taxon_canonical_name=Martinelli@idData$taxon_canonical_name
Martinelli@data$individual_local_identifier=Martinelli@idData$local_identifier
Martinelli <- data.frame(Martinelli@data)

Mimi <- getMovebankData(study=1120749252 , animalName="Mimi 4660", login=loginStored)
Mimi@data$individual_taxon_canonical_name=Mimi@idData$taxon_canonical_name
Mimi@data$individual_local_identifier=Mimi@idData$local_identifier
Mimi <- data.frame(Mimi@data)

Molly <- getMovebankData(study=1120749252 , animalName="Molly 5770", login=loginStored)
Molly@data$individual_taxon_canonical_name=Molly@idData$taxon_canonical_name
Molly@data$individual_local_identifier=Molly@idData$local_identifier
Molly <- data.frame(Molly@data)

Norah <- getMovebankData(study=1120749252 , animalName="Norah 4655", login=loginStored)
Norah@data$individual_taxon_canonical_name=Norah@idData$taxon_canonical_name
Norah@data$individual_local_identifier=Norah@idData$local_identifier
Norah <- data.frame(Norah@data)

Olga <- getMovebankData(study=1120749252 , animalName="Olga 4657", login=loginStored)
Olga@data$individual_taxon_canonical_name=Olga@idData$taxon_canonical_name
Olga@data$individual_local_identifier=Olga@idData$local_identifier
Olga <- data.frame(Olga@data)

Ornette <- getMovebankData(study=1120749252 , animalName="Ornette 4669", login=loginStored)
Ornette@data$individual_taxon_canonical_name=Ornette@idData$taxon_canonical_name
Ornette@data$individual_local_identifier=Ornette@idData$local_identifier
Ornette <- data.frame(Ornette@data)

Peter_Nelson <- getMovebankData(study=1120749252 , animalName="Peter Nelson 5774", login=loginStored)
Peter_Nelson@data$individual_taxon_canonical_name=Peter_Nelson@idData$taxon_canonical_name
Peter_Nelson@data$individual_local_identifier=Peter_Nelson@idData$local_identifier
Peter_Nelson <- data.frame(Peter_Nelson@data)

Pliny <- getMovebankData(study=1120749252 , animalName="Pliny 4675", login=loginStored)
Pliny@data$individual_taxon_canonical_name=Pliny@idData$taxon_canonical_name
Pliny@data$individual_local_identifier=Pliny@idData$local_identifier
Pliny <- data.frame(Pliny@data)

Ripley <- getMovebankData(study=1120749252 , animalName="Ripley 4650", login=loginStored)
Ripley@data$individual_taxon_canonical_name=Ripley@idData$taxon_canonical_name
Ripley@data$individual_local_identifier=Ripley@idData$local_identifier
Ripley <- data.frame(Ripley@data)

Ripley2 <- getMovebankData(study=1120749252 , animalName="Ripley 5771", login=loginStored)
Ripley2@data$individual_taxon_canonical_name=Ripley2@idData$taxon_canonical_name
Ripley2@data$individual_local_identifier=Ripley2@idData$local_identifier
Ripley2 <- data.frame(Ripley2@data)

Riwaka <- getMovebankData(study=1120749252 , animalName="Riwaka 4669", login=loginStored)
Riwaka@data$individual_taxon_canonical_name=Riwaka@idData$taxon_canonical_name
Riwaka@data$individual_local_identifier=Riwaka@idData$local_identifier
Riwaka <- data.frame(Riwaka@data)

Sahti <- getMovebankData(study=1120749252 , animalName="Sahti 4693", login=loginStored)
Sahti@data$individual_taxon_canonical_name=Sahti@idData$taxon_canonical_name
Sahti@data$individual_local_identifier=Sahti@idData$local_identifier
Sahti <- data.frame(Sahti@data)

Sofie <- getMovebankData(study=1120749252 , animalName="Sofie 4674", login=loginStored)
Sofie@data$individual_taxon_canonical_name=Sofie@idData$taxon_canonical_name
Sofie@data$individual_local_identifier=Sofie@idData$local_identifier
Sofie <- data.frame(Sofie@data)

Thelonious <- getMovebankData(study=1120749252 , animalName="Thelonious 4668", login=loginStored)
Thelonious@data$individual_taxon_canonical_name=Thelonious@idData$taxon_canonical_name
Thelonious@data$individual_local_identifier=Thelonious@idData$local_identifier
Thelonious <- data.frame(Thelonious@data)

TonyStark <- getMovebankData(study=1120749252 , animalName="Tony Stark 4659", login=loginStored)
TonyStark@data$individual_taxon_canonical_name=TonyStark@idData$taxon_canonical_name
TonyStark@data$individual_local_identifier=TonyStark@idData$local_identifier
TonyStark <- data.frame(TonyStark@data)

Valoy <- getMovebankData(study=1120749252 , animalName="Valoy 5766", login=loginStored)
Valoy@data$individual_taxon_canonical_name=Valoy@idData$taxon_canonical_name
Valoy@data$individual_local_identifier=Valoy@idData$local_identifier
Valoy <- data.frame(Valoy@data)

Veruca <- getMovebankData(study=1120749252 , animalName="Veruca 4690", login=loginStored)
Veruca@data$individual_taxon_canonical_name=Veruca@idData$taxon_canonical_name
Veruca@data$individual_local_identifier=Veruca@idData$local_identifier
Veruca <- data.frame(Veruca@data)

Vielle <- getMovebankData(study=1120749252 , animalName="Vielle 4670", login=loginStored)
Vielle@data$individual_taxon_canonical_name=Vielle@idData$taxon_canonical_name
Vielle@data$individual_local_identifier=Vielle@idData$local_identifier
Vielle <- data.frame(Vielle@data)

Zola <- getMovebankData(study=1120749252 , animalName="Zola 5212", login=loginStored)
Zola@data$individual_taxon_canonical_name=Zola@idData$taxon_canonical_name
Zola@data$individual_local_identifier=Zola@idData$local_identifier
Zola <- data.frame(Zola@data)

###Consolidate data into one dataframe
data <- rbind(Abby1,Abby2,Avery,BenBob,
              Bob,Bonnie,Carlsberg, Chibi,
              Chloe,Clementina,Da_Vinci,         
              Eli,Ellie,Emma,Fonta_Flora,
              Galena,Gamer,Gillian,Golliath,
              Greg,Ibeth,Inez,Jeff,
              Judy,Kyle,Limon,Mario,
              Martinelli,Mimi,Molly,Norah,
              Olga,Ornette,Peter_Nelson,
              Pliny,Ripley,Ripley2,Riwaka,
              Sahti,Sofie,Thelonious,TonyStark,
              Valoy,Veruca,Vielle,Zola)


###reconstruct data to continuous time
###Model fitting can be slow, may consider running on cluster in future
library(lubridate)
library(ctmm)
data2=as.telemetry(data) 
dt=c(240) #Set sampling rate for variogram
control <- list(method="pNewton",cores=-1)
PROTO <- ctmm(error=FALSE,circle=FALSE)


variograms=lapply(1:50,function(i) variogram(data2[[i]],dt=dt)) #calculate variogram for each animal 
GUESS=lapply(1:50,function(i) ctmm.guess(data2[[i]], CTMM=PROTO,variogram=variograms[[i]],interactive=FALSE)) ##Guess initial parameter values for model fitting
FITS=lapply(1:50,function(i) {
  print(i)
  ctmm.select(data2[[i]],GUESS[[i]],control=control)}) ##Fit best model
library(data.table)

###Simulate from best fit model to interpolate unsampled times
###Average accross 10 simulations for final interpolated track
###Simulation can be slow, may consider running on cluster in future
sims=c()
SIMS=c()
for(i in 1:46){
  sims=c()
  for(j in 1:10){
    print(paste("individual",i,"sim",j,sep="_"))
    SIM <- simulate(data=data2[[i]],object=FITS[[i]], res=1,complete=TRUE,dt=1) #Simulate from movement model at 1 second intervals
    sim=data.frame(SIM)
    index1=which(names(sim)=="longitude")
    index2=which(names(sim)=="latitude")
    sim <- SpatialPointsDataFrame(coords = sim[,c(index1,index2)], data = sim,
                                  proj4string=CRS("+proj=longlat +datum=WGS84"))
    sim <- spTransform(sim, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
    sim=data.frame(sim)
    sims[j]=list(sim)
  }
  sims10=do.call(rbind,sims[c(1:10)]) ##Pull all the simulations for an individual together
  sims10=as.data.frame(sims10)
  sims10=sims10[order(sims10$timestamp),] ##sorts by timestamp
  DT <- data.table(sims10)
  DT$z=DT$longitude.1+1i*DT$latitude.1 ## make coordinates complex numbers for convenience 
  averagesim10=DT[, mean(z,na.rm=T), by = timestamp] #Average all simulations for an individual, summarizes data with mean of each timestamp
  averagesim10=as.data.frame(averagesim10)
  averagesim10$ID=data2[[i]]@info$identity ##annotates all locations with the indiviual being analyzed, for later refenerence when individual are combined
  averagesim10$X=Re(averagesim10$V1) #gets x (real) component of complex mean of locations
  averagesim10$Y=Im(averagesim10$V1) #gets y (imaginary) component of complex mean of locations
  SIMS[i]=list(averagesim10)
}

###calculate grid resolutions for each individual
resolutions=c()

for(i in 1:46){
  print(i)
  if(length(FITS[[i]]$tau)<2){ #in the case where the fitted movement model does not include correllated velocities
    speed=mean(Mod(diff(SIMS[[i]]$V1)),na.rm=T)*60 #Speed set to mean steplength
    Tau_v=1
  }
  else{
    Tau_v=FITS[[1]]$tau[[2]] #When velocities are correlated
    speed=summary(FITS[[i]], units=FALSE)[[3]][which(grepl("speed", rownames(summary(FITS[[i]])[[3]]), fixed=TRUE)),2]
  }
  res=Tau_v*speed #The size of a grid cell is equal to the average speed of the stationary movement process multiplied by the time it takes for the movement process to be uncorrelated, ie the distance between independently used spaces
  resolutions[i]=res
}
detach("package:ctmm", unload = TRUE)

######################################create grid cells#######################################################
library(sp)
library(raster)
library(rgdal)


Gridded=c()
count=1
for(q in 1:50){
  
  route_track <- SIMS[[q]]
  route_track$orientations = c(NA,Arg(diff(route_track$V1))) ##Calculates orientations (head directions)
  route_track$steps = c(NA,Mod(diff(route_track$V1))) ## calculates step-lengths 
  
  colnames(route_track) <- c("timestamp","Z","ID", "X", "Y","orientations", "steps")
  
  route_track_to_grid=SpatialPointsDataFrame(coords=route_track[,c(4,5)],data=route_track)
  #Figure out the grid extents
  e2 <- extent(route_track_to_grid)
  e2@xmin=floor(e2@xmin-resolutions[q])
  e2@xmax=ceiling(e2@xmax+resolutions[q])
  e2@ymin=floor(e2@ymin-resolutions[q])
  e2@ymax=ceiling(e2@ymax+resolutions[q])
  ##########################
  
  route_track_to_grid2 <- raster(e2,res=resolutions[q])
  route_track_to_grid2[] <- 1:ncell(route_track_to_grid2)
  
  
  pol <- rasterToPolygons(route_track_to_grid2) ##Turns the grid into a spatial object
  pols[q]=list(pol) #store in list for later to easily call the grid of a particular individual
  all_dat <- route_track_to_grid %over% pol ###all cell numbers
  route_track$cell=all_dat$layer
  Gridded[q]=list(route_track)
}
Gridded2=do.call(rbind,Gridded)
Gridded2=data.frame(Gridded2)
##Export for DET calculation
write.csv(Gridded2,file="//10.126.19.90/EAS_ind/salavi/simulations/route_data_subset_gridded.csv")


###################################################################ID routes######################################

##Point to DET calculation files
DETs=list.files("C:/Users/salavi/Documents/DETs_by_Cell_Route_Sims3",pattern=".csv", full.names=TRUE)


require(plyr)
require(dplyr)
require(tidyr)
library(ggplot2)
library(ggpubr)
library(mixR)

count=1
calculations=c()
###function to estimate number of modes in distribution, adapted from stackoverflow
get.modes2 <- function(x,adjust,signifi,from,to) {  
  den <- density(x, kernel=c("gaussian"),adjust=adjust,from=from,to=to)
  den.s <- smooth.spline(den$x, den$y, all.knots=TRUE, spar=0.1)
  s.1 <- predict(den.s, den.s$x, deriv=1)
  s.0 <- predict(den.s, den.s$x, deriv=0)
  den.sign <- sign(s.1$y)
  a<-c(1,1+which(diff(den.sign)!=0))
  b<-rle(den.sign)$values
  df<-data.frame(a,b)
  df = df[which(df$b %in% -1),]
  modes<-s.1$x[df$a]
  density<-s.0$y[df$a]
  df2<-data.frame(modes,density)
  df2$sig<-signif(df2$density,signifi)
  df2<-df2[with(df2, order(-sig)), ] 
  #print(df2)
  df<-as.data.frame(df2 %>% 
                      mutate(m = min_rank(desc(sig)) ) %>% #, count = sum(n)) %>% 
                      group_by(m) %>% 
                      summarize(a = paste(format(round(modes,2),nsmall=2), collapse = ',')) %>%
                      spread(m, a, sep = ''))
  colnames(df)<-paste0("m",1:length(colnames(df)))
  return(df)
}


for(q in 1:46){
  print(paste("Animal", q,sep="_"))
  cells=as.character(unique(Gridded[[q]]$cell))
  print(paste("calculate summary stistics for cells of Animal", q,sep="_"))
  
  disresults=c()
  numpts=c()
  nummods=c()
  adjnums=c()
  vars=c()
  ranges=c()
  meanstep=c()
  SDs=c()
  sdstep=c()
  altdist=runif(100000000,min=-pi,max=pi)
  
  #Start progress bar#
  pb <- winProgressBar(title="Calculate cell stistics", label="0% done", min=0, max=length(cells), initial=0)
  
  for(i in 1:length(cells)){
    values=as.numeric(na.omit(Gridded[[q]]$orientations[which(as.character(Gridded[[q]]$cell)==cells[i])]))
    steps=as.numeric(na.omit(Gridded[[q]]$steps[which(as.character(Gridded[[q]]$cell)==cells[i])]))
    vars[i]=var(values,na.rm=TRUE)
    ranges[i]=max(values,na.rm=TRUE)-min(values,na.rm=TRUE)
    meanstep[i]=mean(steps,na.rm=TRUE)
    sdstep[i]=sd(steps)
    
    #calculate hellinger distance
    if(length(values)<2){
      hellinger_tst=NA
    }else{
      hellinger_tst=try(statip::hellinger(values,altdist,-Inf, Inf,method = 2))
    }

    if(length(values)<5){
      modes=0
      mode_count=0
    } else{
      adjnum=0.5
      #calculate number of modes in distribution of orientations
      modes=get.modes2(values,adjust=adjnum,7,min(values)-1,max(values)+1) ##Function can overfit and wildly overestimate modes
      ##If overfit, use mixture model to estimate number of modes
      mode_count=length(modes)
      values2<-data.frame(values=values)
         if(length(modes)>10){
        invisible(capture.output({
          ##Number of modes estimated by using model selection
          ##Fits models of mixtures of up to 10 distributions
          ##Number of modes estimated based on the best fit model as per BIC
          ret1 <- try(select(values, ncomp = 1:10,family = c("normal"),init.method = c("hclust")));
        }))
           mode_count=ret1$ncomp[which(ret1$best=="*")]
      }
     
    }

    numpts[i]=length(values)
    if(length(values)<5){nummods[i]=0}else{
      nummods[i]=mode_count
    }

    adjnums[i]=adjnum
    if(length(values)<2){SDs[i]=0}else{
      SDs[i]=sd(values,na.rm=TRUE)
      }
    res=c(cells[i],numpts[i],hellinger_tst,nummods[i],SDs[i],ranges[i],meanstep[i],sdstep[i])
    disresults[i]=list(res)
    info <- sprintf("%d%% done", round((i/length(cells))*100))
    setWinProgressBar(pb, i/(length(cells))*100, label=info)
    
  }
  close(pb)
  
  disresults2 <- do.call(rbind, disresults)
  disresults2=as.data.frame(disresults2)
  colnames(disresults2)=c("Cell","Density","H_statistic","Modes","sd","spread","mean_steplength","SD_steplength")
  
  print(paste("calculate cell similarity for Animal", q,sep="_"))
  disresults2$Density=as.numeric(as.character(disresults2$Density))
  disresults2$H_statistic=as.numeric(as.character(disresults2$H_statistic))
  disresults2$Modes=as.numeric(as.character(disresults2$Modes))
  disresults2$sd=as.numeric(as.character(disresults2$sd))
  disresults2$spread=as.numeric(as.character(disresults2$spread))
  disresults2$mean_steplength=as.numeric(as.character(disresults2$mean_steplength))
  disresults2$SD_steplength=as.character(disresults2$SD_steplength)
  disresults2$Cell <- as.character(disresults2$Cell)
  
  require(RANN)
  centroids=coordinates(pols[[q]])
  closesttest=nn2(centroids,centroids, k=9)
  closesttest=closesttest$nn.idx
  closesttest=as.data.frame(closesttest)
  ind=c()
  focal=c()
  simil=c()
  
  disresults2$Density_simil <- NA
  disresults2$H_statistic_simil <- NA
  disresults2$Modes_simil <- NA
  disresults2$spread_simil <- NA
  disresults2$sd_simil <- NA
  disresults2$mean_steplength_simil <- NA
  disresults2$SD_steplength_simil <- NA
  disresults2$missing_neighbors <- NA
  pb <- winProgressBar(title="Calculate cell similarity", label="0% done", min=0, max=nrow(disresults2), initial=0)

  ###Calculate the similarity of each cell to it's neighbors (mean squared difference)
  for(i in 1:nrow(disresults2)){
    ind=which(as.character(closesttest[,1])==as.character(disresults2$Cell[i]))
    focal=disresults2[which(disresults2$Cell==as.character(disresults2$Cell[i])),]
    neighbor1=disresults2[which(disresults2$Cell==as.character(closesttest[ind,2])),]
    neighbor2=disresults2[which(disresults2$Cell==as.character(closesttest[ind,3])),]
    neighbor3=disresults2[which(disresults2$Cell==as.character(closesttest[ind,4])),]
    neighbor4=disresults2[which(disresults2$Cell==as.character(closesttest[ind,5])),]
    neighbor5=disresults2[which(disresults2$Cell==as.character(closesttest[ind,6])),]
    neighbor6=disresults2[which(disresults2$Cell==as.character(closesttest[ind,7])),]
    neighbor7=disresults2[which(disresults2$Cell==as.character(closesttest[ind,8])),]
    neighbor8=disresults2[which(disresults2$Cell==as.character(closesttest[ind,9])),]
    lengths=c(length(which(disresults2$Cell==as.character(closesttest[ind,2]))),length(which(disresults2$Cell==as.character(closesttest[ind,3]))),length(which(disresults2$Cell==as.character(closesttest[ind,4]))),length(which(disresults2$Cell==as.character(closesttest[ind,5]))),length(which(disresults2$Cell==as.character(closesttest[ind,6]))),length(which(disresults2$Cell==as.character(closesttest[ind,7]))),length(which(disresults2$Cell==as.character(closesttest[ind,8]))),length(which(disresults2$Cell==as.character(closesttest[ind,9]))))
    missing_neighbors=length(which(lengths==0))
    if(lengths[1]==0){neighbor1[1,]<-rep(0,length(neighbor1))}else{}
    if(lengths[2]==0){neighbor2[1,]<-rep(0,length(neighbor2))}else{}
    if(lengths[3]==0){neighbor3[1,]<-rep(0,length(neighbor3))}else{}
    if(lengths[4]==0){neighbor4[1,]<-rep(0,length(neighbor4))}else{}
    if(lengths[5]==0){neighbor5[1,]<-rep(0,length(neighbor5))}else{}
    if(lengths[6]==0){neighbor6[1,]<-rep(0,length(neighbor6))}else{}
    if(lengths[7]==0){neighbor7[1,]<-rep(0,length(neighbor7))}else{}
    if(lengths[8]==0){neighbor8[1,]<-rep(0,length(neighbor8))}else{}
    
    neighbors=rbind(neighbor1,neighbor2,neighbor3,neighbor4,neighbor5,neighbor6,neighbor7,neighbor8)
    focal2=focal
    disresults2$Density_simil[i]<-mean((as.numeric(na.omit(as.numeric(focal2$Density)-as.numeric(neighbors$Density))))^2,na.rm=TRUE)
    disresults2$H_statistic_simil[i]<-mean((as.numeric(na.omit(as.numeric(focal2$H_statistic)-as.numeric(neighbors$H_statistic))))^2,na.rm=TRUE)
    disresults2$Modes_simil[i]<-mean((as.numeric(na.omit(as.numeric(focal2$Modes)-as.numeric(neighbors$Modes))))^2,na.rm=TRUE)
    disresults2$spread_simil[i]<-mean((as.numeric(na.omit(as.numeric(focal2$spread)-as.numeric(neighbors$spread))))^2,na.rm=TRUE)
    disresults2$sd_simil[i]<-mean((as.numeric(na.omit(as.numeric(focal2$sd)-as.numeric(neighbors$sd))))^2,na.rm=TRUE)
    disresults2$mean_steplength_simil[i]<-mean((as.numeric(na.omit(as.numeric(focal2$SD_steplength)-as.numeric(neighbors$SD_steplength))))^2,na.rm=TRUE)
    disresults2$SD_steplength_simil[i]<-mean((as.numeric(na.omit(as.numeric(focal2$mean_steplength)-as.numeric(neighbors$mean_steplength))))^2,na.rm=TRUE)
    disresults2$missing_neighbors[i] <- missing_neighbors
    setWinProgressBar(pb, i/(nrow(disresults2))*100, label=info)
    
  }
  close(pb)
  
  #########Stop here if you haven't loaded DET data yet#############
  DET=read.csv(DETs[q]) ##Read in DET calculations and add them to data frame
  disresults2=merge(disresults2, DET, by.x = c("Cell"), by.y = c("cell"), all.x = TRUE)
  disresults2$ID=unique(Gridded[[q]]$ID)
  calculations[q]=list(disresults2)
  print(paste("Completed cell calculations for Animal", q,sep="_"))
}


calculations2=do.call(rbind,calculations)
disresults3=calculations2
disresults3$recursions[which(is.na(disresults3$recursions))]=0


for(i in 1:length(cols)){
  disresults3=disresults3[complete.cases(disresults3[,cols[i]]),]
}
##Log density and determinism metrics to make them easier to work with
disresults3$spread_simil[which((is.infinite(disresults3$spread_simil)))]=0
disresults3$log_Density=log(disresults3$Density+.001)
disresults3$log_H_statistic=log(disresults3$H_statistic+.001)
disresults3$log_Density_simil=log(disresults3$Density_simil+.001)
disresults3$log_H_statistic_simil=log(disresults3$H_statistic_simil+.001)
disresults3$log_visits=log(disresults3$visits+.001)
disresults3$log_recursions=log(disresults3$recursions+.001)
disresults3$log_repeats=log(disresults3$repeats+.001)

##Column numbers to go into mixture model
cols=c(4:6,11:13,16,25:31)

###Clustering!
library(ClusterR)

##Find optimal number of clusters
##Optimal number is when BIC levels off
opt_gmm = Optimal_Clusters_GMM(disresults3[,cols], 20, criterion = "BIC", plot_data = TRUE)

##Fit Gaussian mixture model with optimal number of clusters
gmm = GMM(disresults3[,cols], 10, "eucl_dist", "random_subset", 50, 50)
pr = predict_GMM(disresults3[,cols], gmm$centroids, gmm$covariance_matrices, gmm$weights)
disresults3$gmm_label=pr$cluster_labels
testkcalculations_split=split(disresults3,as.factor(disresults3$ID))

##annotate original tracks with cluster labels 
###Will take a while, in future will parallelize
count=1
for(q in 1:46){ 
  print(unique(Gridded[[q]]$ID))
  print(nrow(testkcalculations_split[[q]]))
  
  for(i in 1:nrow(testkcalculations_split[[q]])){
    Gridded[[q]][which(as.character(Gridded[[q]]$cell)==as.character(testkcalculations_split[[q]]$Cell[i])),ncol(Gridded[[q]])+1]=testkcalculations_split[[q]]$gmm_label[i]

  }
  gc()
  colnames(Gridded[[q]])[ncol(Gridded[[q]])]=c("gmm_label")

}



##Plot important features
densplot=ggplot(disresults3)+geom_boxplot(aes(x=as.factor(gmm_label),y=log(Density),fill=as.factor(gmm_label)))+
  ggtitle("a)")+ylab("Density (log)")+xlab(" ")+ labs(fill='Cluster')+theme_classic()
recplot=ggplot(disresults3)+geom_boxplot(aes(x=as.factor(gmm_label),y=log(recursions+.001),fill=as.factor(gmm_label)))+
  ggtitle("b)")+ylab("Recursions (log)")+xlab(" ")+ labs(fill='Cluster')+theme_classic()

repplot=ggplot(disresults3)+geom_boxplot(aes(x=as.factor(gmm_label),y=log(repeats+.001),fill=as.factor(gmm_label)))+
  ggtitle("c)")+ylab("Repeats (log)")+xlab(" ")+ labs(fill='Cluster')+theme_classic()

ggarrange(densplot,recplot,repplot,ncol=3,common.legend = TRUE, legend = "bottom")


###Calculate routineness score

meansden=with(disresults3, tapply(log_Density, as.factor(gmm_label), mean)) #mean of each cluster 
meansdenQuant=cut(meansden, ##Break into quartiles
                  breaks = qu <- quantile(meansden),
                  labels = names(qu)[-1], include.lowest=TRUE)
denquant=cbind(names(meansden),as.character(meansdenQuant))
meanrep=with(disresults3, tapply(log_repeats, as.factor(gmm_label), mean))#mean of each cluster 
meansrepQuant=cut(meanrep,##Break into quartiles
                  breaks = qu <- quantile(meanrep),
                  labels = names(qu)[-1], include.lowest=TRUE)
repquant=cbind(names(meanrep),as.character(meansrepQuant))

meanrecur=with(disresults3, tapply(log_recursions, as.factor(gmm_label), mean))#mean of each cluster 
meansrecurQuant=cut(meanrecur,##Break into quartiles
                    breaks = qu <- quantile(meanrecur),
                    labels = names(qu)[-1], include.lowest=TRUE)
recurquant=cbind(names(meanrecur),as.character(meansrecurQuant))

denquant=data.frame(denquant)
colnames(denquant)=c("Cluster","denquant")
denquant$denquant=as.numeric(sub("%", "", denquant$denquant))/100 #reformat into decimal

repquant=data.frame(repquant)
colnames(repquant)=c("Cluster","repquant")
repquant$repquant=as.numeric(sub("%", "", repquant$repquant))/100 #reformat into decimal

recurquant=data.frame(recurquant)
colnames(recurquant)=c("Cluster","recurquant")
recurquant$recurquant=as.numeric(sub("%", "", recurquant$recurquant))/100 #reformat into decimal

denreq=merge(denquant,repquant, by="Cluster")
denreq=merge(denreq,recurquant, by="Cluster")

denreq$Label=(denreq$denquant*denreq$repquant*denreq$recurquant)/sum(denreq$denquant*denreq$repquant*denreq$recurquant)

denreq$Cluster=as.numeric(denreq$Cluster)                

#Add routineness score of each cluster to data frame
for(j in 1:10){
  disresults3$Label[which(disresults3$gmm_label==denreq$Cluster[j])]=denreq$Label[(which(denreq$Cluster==j))]
  
}


###########Prepare data for plotting####################
###Plot only real fixes, not interpolated

  data3 <- SpatialPointsDataFrame(coords = data[,c(20,19)], data = data,
                                  proj4string=CRS("+proj=longlat +datum=WGS84"))
  data3 <- spTransform(data3, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
  data3=data.frame(data3)
  data3$individual_local_identifier=as.factor(data3$individual_local_identifier)
  splitdata=split(data3,data3$individual_local_identifier)

  finaldata=c()
  for(i in 1:46){
    splitdata[[i]]$individual_local_identifier=as.character(splitdata[[i]]$individual_local_identifier)
    splitdata[[i]]$X2= splitdata[[i]]$location_long.1+1i* splitdata[[i]]$location_lat.1
    temp=merge(splitdata[[i]], Gridded[[i]], by.x = c("timestamp","individual.local.identifier",
                                                      "individual.taxon.canonical.name" ), 
               by.y = c("timestamp","ID","species"), all.x = TRUE)
    temp=temp[complete.cases(temp$gmm_label),]
    finaldata[i]=list(temp)
  }
  
  for(i in 1:46){
    for(j in 1:10){
      #finaldata[[i]]$Label=NA
      finaldata[[i]]$Label[which(finaldata[[i]]$gmm_label==denreq$Cluster[j])]=denreq$Label[(which(denreq$Cluster==j))]
      
    }
  }
  plotdata=finaldata[c(9,4,11,12)] #Pick one of each species to plot
  plotdata=do.call(rbind,plotdata)
  plotdata$individual_local_identifier=as.factor(plotdata$individual_local_identifier)
  
  ##Plot routineness
  ggplot(plotdata[complete.cases(plotdata$gmm_label),])+geom_path(aes(X.y,Y,color=Label, group = 1))+scale_color_gradientn(colours = colorspace::heat_hcl(7))+theme_classic()+
    theme(legend.position="bottom") + facet_wrap(~individual_taxon_canonical_name, ncol=2, scales = "free")+ theme(aspect.ratio = 1)+ylab("Northing")+xlab("Easting")+
    labs(color="Routineness score")
  
  ##Plot tracks on island
  ggplot()+geom_polygon(data=outline_df,aes(x=long,y=lat,group=group),color="black",fill=NA)+
    geom_path(data=plotdata[complete.cases(plotdata$gmm_label),],aes(X.y,Y, group=individual_taxon_canonical_name, color=individual_taxon_canonical_name, alpha=1))+ coord_fixed()+
    theme_classic()+ theme(legend.position="bottom",legend.title = element_blank()) + xlab("Easting") + ylab("Northing")+annotation_scale()+ scale_alpha(guide = 'none')
  
  ##Plot cell examples
  par(mfrow=c(1,3))
  plot(pols[[12]][221,],col="red", main="a)")
  lines(Gridded[[12]]$X2)
  plot(pols[[12]][111,],col="red", main="b)")
  lines(Gridded[[12]]$X2)
  plot(pols[[12]][85,],col="red", main="c)")
  lines(Gridded[[12]]$X2)
  

#####Compare routineness of species 
library(rstan)
library(brms)
library(ggplot2)
options(mc.cores = parallel::detectCores()) ##parallelize the chains

model1=brm(Label~Species+(1|ID),data=disresults3, chains=2)
summary(model1)
plot(model1)
viz=conditional_effects(model1,spaghetti=FALSE)
viz=plot(viz, plot=FALSE)
viz=viz[[1]]+ylab("Routineness score")+ggtitle("a)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_classic()
viz=viz+theme(axis.text.x = element_text(angle = 45, hjust = 1))

results=summary(model1)
results=data.frame(results$fixed)
results$covariate=c("Intercept", "Cebus capucinus", "Nasua narica", "Potos flavus")
posteriors=ggplot(results,aes(x=Estimate,y=covariate))+geom_point()+geom_linerange(aes(xmin=l.95..CI,xmax=u.95..CI))+
  geom_vline(xintercept = 0)+
  theme_classic() +ylab("Covariate")+ggtitle("b)")
posteriors
ggarrange(viz,posteriors)
