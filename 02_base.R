
# grids sampled ------------ 

# sample sizes
phase_nGrids<- rowSums(table(ph@other$ind.metrics$phaseNo, ph@other$ind.metrics$gridId) > 0)
phase_nGrids
# grids 
    
grids_more2n <- names(which(colSums(table(ph@other$ind.metrics$phaseNo, 
                                        ph@other$ind.metrics$gridId)) > 2))
grids_2n <- names(which(colSums(table(ph@other$ind.metrics$phaseNo, 
                                         ph@other$ind.metrics$gridId)) == 2))
gridNames <- c(grids_more2n, grids_2n[3])
gridNames %>% length # the grids to simulate and to be used in the analysis 

# base of simulation ---------

ph2 <- gl.keep.pop(ph, as.pop = 'phaseNo', pop.list = 'I1')

nSample <- table(ph@other$ind.metrics$phaseNo,
                 factor(ph@other$ind.metrics$gridId)) %>% 
  data.frame()

gridsKeep <- nSample$Var2[nSample$Var1 == 'I1' & nSample$Freq>3] %>% unique()
gridsKeep %>% length # 24 alf grids

pop(ph2) <- ph2@other$ind.metrics$gridId
glBase <- gl.keep.pop(ph2, pop.list = gridsKeep)
glBase <- gl.filter.monomorphs(glBase)

gl.report.monomorphs(glBase)

popsReal <- seppop(glBase); popsReal <- popsReal[c(1:24,1:24)]

names(popsReal)[25:48] <- gridNames[!gridNames %in% gridsKeep] 

for(i in 1:length(popsReal)){
  popsReal[[i]]@other$ind.metrics$gridId <- names(popsReal)[i]
  pop(popsReal[[i]]) <- popsReal[[i]]@other$ind.metrics$gridId
}

sapply(popsReal, function(x) table(pop(x)))
# step 2: define site gene flow probabilities ----------

gridNames %>% length
xy <- matrix(nrow = length(gridNames),
             ncol = 2,
             data = 1:(length(gridNames)*2))
rownames(xy) <- gridNames

mat <- as.matrix(dist(xy))
mat[mat != 0] <- 1/(nrow(mat)-1) # equal dispersal

colSums(mat)
nrow(mat)
mat[1:5, 1:5]
dim(mat)
# step 3: simulate alf ---------
alfreal <- gl.alf(glBase)

popSim <- mclapply(popsReal, function(x) gl.sim.ind(x, n = 100, x@pop[1]), mc.cores = 24)

popSim <- lapply(popSim, fxsex)
simStart <- reduce(popSim, em.gl.join)
pop(simStart) %>% table


# compare alf --------------------
alfsim <- gl.alf(simStart)

alfsim$type <- 'sim'
alfreal$type <- 'real'

alfdiff <- alfreal$alf1 - alfsim$alf1 
hist(alfdiff)
max(abs(alfdiff))

rbind(alfreal, alfsim) %>% 
  ggplot(aes(x = alf2, fill = type))+
  geom_histogram(bins = 30, position = position_dodge())+
  theme_classic()
