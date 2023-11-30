
# gen zero ---------------------------------------------------------------------
em.gl.indmetrics <- function(glsim, gen = 0, disp.rate = 0) {
  genn <- paste0("gen", sprintf("%02d", gen))
  glsim@other$ind.metrics <- data.frame(id = glsim@ind.names, 
                                        pop = glsim@pop,
                                        sex = glsim@other$sex,
                                        generation = 0)
  return(glsim@other$ind.metrics)
}  

# add sex if missing
fxsex <- function(x){
  x@other$sex <- sample(c('male', 'female'), nInd(x), replace = T)
  
  x@other$ind.metrics <- em.gl.indmetrics(x)
  x@other$ind.metrics$id <- paste0(x@other$ind.metrics$pop, '_', 
                                   x@other$ind.metrics$id)
  return(x)
}

# join two genlight objects and their ind.metric data
em.gl.join <- function(x1, x2){
  x <- rbind(x1, x2)
  x@other$ind.metrics <- rbind(x1@other$ind.metrics, x2@other$ind.metrics)
  return(x)
}

# simulate generation zero 
gen.zero.offspring <- function(pstart, site, n.offspring = 4){
  dads <- gl.keep.pop(pstart, as.pop = "sex", pop.list = "m")
  mums <- gl.keep.pop(pstart, as.pop = "sex", pop.list = "f")
  p <- gl.sim.offspring(fathers = dads, mothers = mums, 
                        noffpermother = n.offspring, sexratio = 0.5)
  n.ind <- nInd(p)
  p@pop <- factor(rep(site, n.ind))
  indNames(p) <- paste0(p@pop, "_", 1:n.ind)
  p@other$ind.metrics <- em.gl.indmetrics(p)
  return(p)
}

# simulate generation zero 
gen.zero.alf <- function(pstart, site, Nsize = 100 ){
  p <- gl.sim.ind(pstart, n = Nsize, popname = site)
  n.ind <- nInd(p)
  p@pop <- factor(rep(site, n.ind))
  indNames(p) <- paste0(p@pop, "_", 1:n.ind)
  p@other$sex <- sample(c('male', 'female'), n.ind, replace=T)
  p@other$ind.metrics <- em.gl.indmetrics(p)
  return(p)
}


em.gl.sim.ind <- function (x, n = 250, popname = NULL) 
{
  p <- as.matrix(x)
  nind <- nInd(x)
  nx <- round(n/nind) # could change this where nx is what you use, eg nx = 10;
                      #                                          n=4>40 n=10>100
  
  simind <- p
  for(i in 1:(nx-1)){
  new <- apply(p, MARGIN = 2, function(x) sample(x, nind, replace = F))
  rownames(new) <- paste0(rownames(p),'_sim', i)
  simind <- rbind(simind, new)
  }
  glsim <- new("genlight", gen = simind, ploidy = 2, 
               ind.names = rownames(simind), loc.names = locNames(x), loc.all = x@loc.all, 
               position = position(x), pop = rep(popname, nrow(simind)))
  return(glsim)
}

# gene flow --------------------------------------------------------------------

gene.flow <- function(gl, probMatrix, mRate){

  zero <- gl
  zero@other$ind.metrics$disperse <- FALSE
  popId <- zero@other$ind.metrics[, c("pop", "id")]
  popId$index <- 1:nrow(popId)
  popSplit <- split(popId, popId$pop)
  
  em.mig <- function(x, rate) {
    Ne <- length(x$index)
    Nm <- Ne*rate
    left <- Nm %% 1
    if(left != 0) Nm <- sample(c(floor(Nm),ceiling(Nm)), 1, prob = c(1-left, left))
    migrants <- sample(x$index, size = Nm)
    # to stop loss of pop (f-m):
    if(Ne < 5) migrants <- sample(x$index, size = 0) 
     return(migrants)
    }

  disperID <- lapply(popSplit, em.mig, rate = mRate) 
  dID <- do.call("c", disperID)
  zero@other$ind.metrics$disperse[dID] <- TRUE
  meta <- zero@other$ind.metrics
  sites <- colnames(probMatrix)
  for(i in 1:nrow(meta)){
    if(meta$disperse[i]){
      si <- meta$pop[i] 
      siDisp <- probMatrix[,si]
      meta$pop[i] <- sample(sites, 1, prob = siDisp)
    }
    
  }
  zero@other$ind.metrics <- meta 
  pop(zero) <- meta$pop
  return(zero)
}

# next genertion ---------------------------------------------------------------
next.gen <- function(popx, n.offspring = 4){
  p <- popx
  site <- unique(p@other$ind.metrics$pop)
  pop(p) <- p@other$ind.metrics$sex
  fmp <- seppop(p)
  nextGen <- gl.sim.offspring(fathers = fmp[["male"]],
                              mothers = fmp[["female"]],
                              noffpermother = n.offspring,
                              sexratio = 0.5)
  prevGen <- p@other$ind.metrics$generation[1]
  meta <- data.frame(id = sub("Po", site, nextGen@ind.names), 
                     pop = site,
                     sex = nextGen@other$sex,
                     disperse = FALSE,
                     generation = prevGen + 1)
  
  pop(nextGen) <- meta$pop
  nextGen@ind.names <- meta$id
  nextGen@other$ind.metrics <- meta
  
  # survival 
 # if(nInd(nextGen)>N.size) nextGen <- nextGen[sample(1:nInd(nextGen), N.size),]
  
  
  return(nextGen)
  
}

the.survivors <- function(nextGen, max.N =500){
  # survival 
   N <- nInd(nextGen)
   if(N>max.N) nextGen <- nextGen[sample(1:N, max.N),]
  
  return(nextGen)
  
}


em.gl.join <- function(x1, x2){
  x <- rbind(x1, x2)
  x@other$ind.metrics <- rbind(x1@other$ind.metrics, x2@other$ind.metrics)
  return(x)
}

# big simulation -------------
em.simulate <- function(zero, mat, conditions){
  sim <- list()
  gen <- zero#sim[[11]]
  
  system.time({
    for(i in 1:nrow(conditions)){
      # setup
      icon <- conditions[i,]
      m.ind <- icon$migration
      max.N <- icon$N 
      pha <- icon$phase
      off <- icon$offspring
      
      # mice dynamics
      gen <- gene.flow(gen, mat, m.ind)
      pops <- seppop(gen)
      genRepro <- lapply(pops, next.gen, n.offspring = off)
      genSurvive <- lapply(genRepro, the.survivors, max.N = max.N)
      gen <- reduce(genSurvive, em.gl.join)
      gen@other$ind.metrics$iteration <- i
      gen@other$ind.metrics$nInd <- nInd(gen)
      gen@other$ind.metrics <- cbind(gen@other$ind.metrics, icon)
      sim[[i]] <- gen
      txt <- paste("  Generation:", i, pha)
      cat(paste0("\033[0;", 36, "m", txt,"\033[0m","\n"))
      
    }
  })
  
  return(sim)
}
