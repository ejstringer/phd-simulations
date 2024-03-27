## run simulation ----------
dir.create('/data/scratch/emily/simulations/sy/')
phaseNo <- paste0(rep(c('L','I', 'D'), 3), rep(1:3, each = 3))[-1]
sim_schedule <- paste0('simulationsy_', sprintf("%02d", 1:20), '.rds')
initialisedSelected <- readRDS('./output/initialised_simulation_gen0_sy.rds')
simdir <- '/data/scratch/emily/simulations/sy/genlights/'
condir <- '/data/scratch/emily/simulations/sy/meta/'
dir.create(simdir)
dir.create(condir)

# no adjustment
mAdjust <- 0
# original adjustment
#mAdjust <- rep(c(0.02,-0.02, 0.02, 0, 0, -0.02, 0, 0), phaselength)

system.time({
for (i in 1:length(sim_schedule)) {
  
  saved_sims <- list.files(path = '/data/scratch/emily/simulations/sy/genlights/')
  if(!sim_schedule[i] %in% saved_sims){
    print(paste('running simulation', i))
    
    # adjust m
    m<- conditions$migration+mAdjust
    
    m <- ifelse(m <0 , 0, m)
    m <- ifelse(m >= 1, 1, m)
    print(round(m,3))
    conditions$migration <- m
    
    
  sim <- em.simulate(initialisedSelected, mat, conditions, ncores = 20)
  
  simJoin <- reduce(sim, em.gl.join)
  simJoin@other$ind.metrics$phaseNo <- factor(simJoin@other$ind.metrics$phaseNo,
                                              levels = phaseNo)
  pop(simJoin) <- simJoin@other$ind.metrics$phaseNo
  
  simPhase <- seppop(simJoin)

  simPhase <- lapply(simPhase, em.gl_assign_pop, define.pop = 'pop')
  
  names(sim) <-  paste0('gen', sprintf("%02d", 1:length(sim)))
  
  simAll <- c(sim, simPhase)
  names(simAll)
    
  
  filename <- sim_schedule[i]
  saveRDS(simAll, paste0(simdir, filename))
  
  conditions_save <- conditions
  conditions_save$nInd <- sapply(sim, nInd)
  conditions_save$perPop <- round(sapply(sim, nInd)/sapply(sim, nPop))
  conditions_save$simulation <- sprintf("%02d", i)
  
  write.csv(conditions_save, paste0(condir, 
                                    sub('\\.rds', '_conditions.csv', filename)),
            row.names = FALSE)
  
  }
}
})




