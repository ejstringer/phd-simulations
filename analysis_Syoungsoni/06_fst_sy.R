# fst

library(parallel)
library(foreach)
library(doParallel)
library(dartR)
# setup dir info --------------------------------------
phaselength <- c(1, 1, 1, 2, 1, 3, 1, 1)
simdir <- '/data/scratch/emily/simulations/sy/genlights/'
fstdir <- '/data/scratch/emily/simulations/sy/fst/'
condir <- '/data/scratch/emily/simulations/sy/meta/'
dir.create(fstdir)

sim_list <- list.files('/data/scratch/emily/simulations/sy/genlights/')


# setup cores -----------------
ncores <- 11+8
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)
tofst <- 1:ncores

# run fst ---------------------
system.time({
for(i in 1:length(sim_list)){
filename <- sim_list[i]
fstfilename <- sub('\\.rds', '_fst.rds', filename)
fst_list <- list.files('/data/scratch/emily/simulations/sy/fst/')
if(!fstfilename %in% fst_list){
  print(paste('calculate fst for simulation', i))
confilename <- paste0(condir, sub('.rds', '_conditions.csv', filename))
simAll <- readRDS(paste0(simdir, filename))

print(paste('finish loading file...', Sys.time()))
  simFst <- foreach(k = tofst) %dopar% {
    x <- TRUE
    while(x){
      fst <- dartR::gl.fst.pop(simAll[[k]], nboots = 1)
      if(is.matrix(fst)) x <- FALSE
    } # end while
    fst
  } # end foreach
  print(paste('finish calculating fsts...', Sys.time()))

saveRDS(simFst, paste0(fstdir, fstfilename))

conditionsx <- read.csv(confilename)
conditionsx$fst <- sapply(simFst[1:11], mean, na.rm = T)
conditionsx$fstphase <- rep(sapply(simFst[12:19], mean, na.rm = T), phaselength)
write.csv(conditionsx, confilename, row.names = FALSE)
}# end if
}# end for
})# end sys time
stopCluster(cl)




