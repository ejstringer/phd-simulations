# all freq ----------
em.alf_df <- function(ph_phase, meta){
  alfphase <- lapply(ph_phase, gl.alf)
  
  for(i in 1:length(alfphase)){
    
    a <- alfphase[[i]]
    alfphase[[i]]$loci <- rownames(a)
    alfphase[[i]]$phaseNo <- names(alfphase)[i]
    
  }
  
  alfdf <- do.call('bind_rows', alfphase) %>% 
    left_join(meta) %>% 
    mutate(npp.log = log(npp),
           uid = sub('-', '_', loci),
           uid = sub('-.*', '', uid),
           uid = sub('_', '-', uid))
  rownames(alfdf) <- NULL
  table(alfdf$phaseNo)
  alfdf$phaseNo <- factor(alfdf$phaseNo,
                          levels = levels(ph@other$ind.metrics$phaseNo))
  
  return(alfdf)
}

# subset -------
samplex <- table(ph_real@other$ind.metrics$phaseNo, 
                 factor(ph_real@other$ind.metrics$gridId))[-1,]
ncol(samplex) 

em.sample_sim <- function(samplex, simxx){
  samplelist <- list()
  
  for(i in 1:nrow(samplex)){
    simxx[[i]]@ind.names <- paste0(simxx[[i]]@ind.names, '_g',
                                   simxx[[i]]@other$ind.metrics$gen)
    metaInd <- simxx[[i]]@other$ind.metrics
    nChosenOnes <- c()
    for (j in 1:ncol(samplex)) {
      n <- samplex[i,j]
      
      if(n > 0){
        p <- rownames(samplex)[i]
        g <- colnames(samplex)[j]
        nChosen <- sample(which(metaInd$phaseNo == p & metaInd$pop == g), n)
        nChosenOnes <- c(nChosenOnes, simxx[[i]]@ind.names[nChosen])
      }
    }
    samplelist[[i]] <- nChosenOnes
  }
  samplex[1:2,1:2]
  rowSums(samplex)
  sapply(samplelist, length)
  
  newSimsub <- list()
  for(x in 1:length(simxx)){
    
    newSimsub[[x]] <- simxx[[x]][which(simxx[[x]]@ind.names %in% samplelist[[x]]),]
    
  }
  newSimsub
  sapply(newSimsub, nInd)
  names(newSimsub) <- names(simxx)
  return(newSimsub)
}


# npp random -----------
meta <- ph@other$ind.metrics %>% 
  group_by(period, phaseNo, phase, npp) %>% 
  summarise(n = n(),
            ngrids = length(unique(gridId))) %>% 
  rename(nppTrue = npp)
random <- meta[,c('phaseNo', 'nppTrue')]  %>% unique()

set.seed(444)
random$npp <- sample(random$nppTrue, 9); cor(random$nppTrue, random$npp)
plot(log(npp)~log(nppTrue), data = random)

meta <- left_join(meta, random)


# real alf -----------


alfReal <- em.alf_df(ph_phase = ph_phase, meta)
dir.create('/data/scratch/emily/simulations/alf')
saveRDS(alfReal, '/data/scratch/emily/simulations/alf/real_alf.rds')

alfdf <- filter(alfReal, phaseNo != 'L1')

# sim alf ------------
# simulation ----


ncores <- 20
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)


sim_list <- list.files('/data/scratch/emily/simulations/genlights/')

alf_nSim <- foreach(i = 1:20) %dopar% {
  library(dartR)
  library(tidyverse)
  
  print(paste('sim', i, '....Go!!'))
  filename <- sim_list[i]
  fileload <- list.files('/data/scratch/emily/simulations/genlights/',
                         full.names = T)[i]
  simx <- readRDS(fileload)[20:27]
  nSim_list <- list()
  
  for(j in 1:50) {
    
    alffilename <- sub('\\.rds', paste0('_alf_n', sprintf("%02d", j), '.rds'),
                       filename)
    alf_list <- list.files('/data/scratch/emily/simulations/alf/')
    if(!alffilename %in% alf_list){
      print(paste('calculate alf for simulation', i,'subsample', j, '....Go!!'))
      
      nsim<- em.sample_sim(samplex, simx)
      system.time(nSimalf <- em.alf_df(nsim, meta))
      nSimalf$i <- i+(j/100)
      saveRDS(nSimalf, paste0('/data/scratch/emily/simulations/alf/', 
                              alffilename))
      
      nSim_list[[j]] <- nSimalf 
    }else{
      nSim_list[[j]] <- readRDS(paste0('/data/scratch/emily/simulations/alf/', 
                                       alffilename))
    }
  }
  
  to_save_nSim <- do.call('bind_rows', nSim_list)
  to_save_nSim
}

stopCluster(cl)

change_nSim2 <- do.call('bind_rows', alf_nSim)
table(change_nSim2$i)

