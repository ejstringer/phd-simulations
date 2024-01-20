## run simulation ----------

phaseNo <- paste0(rep(c('L','I', 'D'), 3), rep(1:3, each = 3))[-1]
sim_schedule <- paste0('simulation_', sprintf("%02d", 1:50), '.rds')
initialisedSelected <- readRDS('./output/initialised_simulation_gen0.rds')
simdir <- '/data/scratch/emily/simulations/genlights/'
condir <- '/data/scratch/emily/simulations/meta/'


system.time({
for (i in 1:length(sim_schedule)) {
  
  saved_sims <- list.files(path = '/data/scratch/emily/simulations/genlights/')
  if(!sim_schedule[i] %in% saved_sims){
    print(paste('running simulation', i))
  sim <- em.simulate(initialisedSelected, mat, conditions, ncores = 30)
  
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




