
# setup --------
source('code/libraries.R')
source('code/functions_simulate_geneflow.R')

# what we know: low pop size = 25* with migration rate 0
#               2007 had an increase rate of ~4.17, round to 4*
#               number of grids to take alf from is 20
#               2007 Fst = 0.014 from Ch1 and 0.017* from Ch2
#               generation time every 3months during booms (increase decrease)
#               geneartion time every year during busts (low)
#               migration = 0 during decreases and lows

# step 1: data -----------
ph <- readRDS('../phd-analysis2/output/Pseudomys_hermannsburgensis_filtered_genotypes.rds')
ph2 <- gl.keep.pop(ph, as.pop = 'trip', pop.list = '2007-09-01')

nSample <- table(ph@other$ind.metrics$trip,
                 factor(ph@other$ind.metrics$gridId))
gridsKeep <- names(nSample[22,][nSample[22,]>3])
gridsKeep %>% length


pop(ph2) <- ph2@other$ind.metrics$gridId
glBase <- gl.keep.pop(ph2, pop.list = gridsKeep)
glBaseMaf <- gl.filter.maf(glBase, threshold = 0.1)

system.time(fstrealmaf <- gl.fst.pop(glBaseMaf, nboots = 1))
mean(fstrealmaf, na.rm = T)
median(fstrealmaf, na.rm = T)

popsReal <- seppop(glBaseMaf)

# step 2: define site gene flow probabilities ----------
xy <- matrix(nrow = length(gridsKeep),
             ncol = 2,
             data = 1:(length(gridsKeep)*2))
rownames(xy) <- gridsKeep

mat <- as.matrix(dist(xy))
mat[mat != 0] <- 1/(nrow(mat)-1) # equal dispersal

colSums(mat)
nrow(mat)
mat

# step 3: simulate alf ---------
popSim <- pblapply(popsReal, function(x) gl.sim.ind(x, n = 100, x@pop[1]))

popSim <- lapply(popSim, fxsex)

popSim$CS2@other$ind.metrics <- lapply(popSim, em.gl.indmetrics)
simStart <- reduce(popSim, em.gl.join)

system.time(fststart <- gl.fst.pop(simStart, nboots = 1))
mean(fststart, na.rm = T)
median(fststart, na.rm = T)

# step 4: initalise simulation ---------

ne <- nInd(popSim$FRN1)
fst <- 0.017
m <- round((1/fst-1)/(4*(ne)), 2)

initialise_conditions <- data.frame(gen = 1:5, 
                                    migration = m, N = ne, 
                                    phase = 'I', offspring = 4,
                                    npop = nPop(simStart)) %>% 
  mutate(leaving = N*migration, Nm = leaving/(npop-1))

initialiseSim <- em.simulate(simStart, mat,initialise_conditions)


# step 5: define simulation conditions ---------




# step 6: run simulation ----------

# step 7: Fsts ---------

