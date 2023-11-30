
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

# system.time(fstrealmaf <- gl.fst.pop(glBaseMaf, nboots = 1))
# mean(fstrealmaf, na.rm = T)
# median(fstrealmaf, na.rm = T)

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
mat[1:5, 1:5]

# step 3: simulate alf ---------
popSim <- pblapply(popsReal, function(x) gl.sim.ind(x, n = 100, x@pop[1]))

popSim <- lapply(popSim, fxsex)


simStart <- reduce(popSim, em.gl.join)

# system.time(fststart <- gl.fst.pop(simStart, nboots = 1))
# mean(fststart, na.rm = T)
# median(fststart, na.rm = T)

# step 4: initalise simulation ---------

ne <- nInd(popSim$FRN1)
fst <- 0.017
mEqu <- ceiling((1/fst-1)/(4*(ne))*100)/100 # at equalibrium
m <- mEqu#*2 

initialise_conditions <- data.frame(gen = 1:10, 
                                    migration = m, N = ne, 
                                    phase = 'I', offspring = 4,
                                    npop = nPop(simStart)) %>% 
  mutate(leaving = N*migration, Nm = leaving/(npop-1))

initialise_conditions$migration[1] <- 0.615

initialiseSim <- em.simulate(simStart, mat,initialise_conditions)

sapply(initialiseSim, nPop)
sapply(initialiseSim, nInd)

fstinit <- pblapply(initialiseSim, gl.fst.pop, nboots = 1)

data.frame(fst = sapply(fstinit, mean, na.rm = T),
           gen = 1:length(fstinit)) %>% 
  ggplot(aes(gen, fst))+
  geom_hline(yintercept = fst, colour = 'grey', lty = 2)+
  geom_point(size =4)+
  theme_classic()

fstinit[[4]] %>% mean(., na.rm = T)
# step 5: define simulation conditions ---------
phaseNo <- paste0(rep(c('I', 'D', 'L'), 3), rep(1:3, each = 3))[-9]

phaselength <- c(2, 2, 2, 4, 2, 4, 2, 2)

conditions <- data.frame(gen = 1:sum(phaselength), period = rep(1:8, phaselength),
                         phaseNo = rep(phaseNo[-1], phaselength), 
                         grids = sample(4:10, sum(phaselength), replace = T)
) %>% 
  mutate(phase = str_sub(phaseNo, 1,1),
         rep = str_sub(phaseNo, 2,2),
         grids = ifelse(phase == 'I', grids + 10, grids),
         grids = ifelse(phase == 'D', grids + 3, grids),
         migration = ifelse(phase == 'I', 0.25, 0),
         N = case_when(
           phase == 'L' ~ 25,
           phase == 'I' ~ 250,
           phase == 'D' ~ 250/2,
         ),
         N = ifelse(phase == 'D' & duplicated(phaseNo), 250/4, N),
         N = ceiling(N),
         offspring = 4) %>% 
  mutate(rep = ifelse(phaseNo == 'L3', 2, rep),
         rep = ifelse(phaseNo == 'L2', 1, rep),
         genx = c(1:6, 1:10, 1:4),
         years = c(1,1.25,2,2.2,3,4,
                   1.1,1.2,1.3,1.4, 2, 2.2, 3,4,5,6,
                   1,1.2,2,2.2))



# step 6: run simulation ----------

# step 7: Fsts ---------

