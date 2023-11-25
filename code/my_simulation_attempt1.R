source('code/libraries.R')
source('code/functions_simulate_geneflow.R')

# genetics ----------
ph2 <- readRDS('../data/simulation_base_phgenlight_2007sep.rds')
pop(ph2) <- ph2@other$ind.metrics$gridId
ph2 <- gl.filter.callrate(ph2, threshold = 1)
ph2 <- gl.filter.monomorphs(ph2)
glBase <- gl.keep.pop(ph2, as.pop = "sex", pop.list = c("m", "f"))
table(glBase@other$ind.metrics$sex, exclude = F)

# spatial ----------
grids <- read.csv("../data/em_gridcoortable.csv") %>% 
  filter(complete.cases(lon))
grids[grids$coor== "site mean",c("lat", "lon")] <- NA
xyGrids <- rgdal::project(as.matrix(grids[,c("lon", "lat")]),
                          proj = '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
rownames(xyGrids) <- grids$gridId
xy <- xyGrids[rownames(xyGrids) %in% c("KSE2", "MC4",  "SS3",  "FRN1", "FRS2", "MC2"),]
#xy <- xyGrids[rownames(xyGrids) %in% c("MC2", "MC4"),]
nrow(xy)

xy[,1] <- rep(1:3, 2)
xy[,2] <- rep(1:2, each = 3)
plot(xy, col = rainbow(6), pch = 16)
lcost <- xy
mat <- as.matrix(dist(lcost))
mat[mat != 0] <- 1/(nrow(mat)-1) # equal dispersal

colSums(mat)
mat
## zero (gen1) ----------

sites <- rownames(xy) # populations
p <- lapply(sites, function(x) gen.zero(glBase, x, 4))
zero <- reduce(p, em.gl.join)
table(pop(zero)) 

fstList <- list()
indkeep <- lapply(1:30, function(x) sample(1:nInd(zero), 6*9))
lapply(indkeep, function(x) zero[x,])
fst <- pbapply::pblapply(indkeep, function(x) gl.fst.pop(zero[x,])$Fsts)
sapply(fst, function(x) log(mean(x[x>0], na.rm = T))) %>% mean %>% exp
zero@other$ind.metrics
n.zero <- zero[indkeep,]
system.time(fst <- gl.fst.pop(n.zero))

lapply(indkeep, )

fst$Fsts[fst$Fsts>0] %>% mean(.,na.rm = T) %>% log

# my system ---------
phaseNo <- paste0(rep(c('L', 'I', 'D'), 3), rep(1:3, each = 3))
100/10
882/231
1076/167
738/353

25*10

phaselength <- c(1, 2, 2, 4, 2, 4, 2, 2)

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
         N = ifelse(phase == 'D' & duplicated(phase), 250/4, N),
         N = ceiling(N))


# sim -----

head(conditions)
mat
sim <- list()
gen <- zero

system.time({
for(i in 1:3){
  # setup
  icon <- conditions[i,]
  m.ind <- icon$migration
  max.N <- icon$N 
  pha <- icon$phase
  
  # mice dynamics
  gen <- gene.flow(gen, mat, m.ind)
  pops <- seppop(gen)
  genRepro <- lapply(pops, next.gen, n.offspring = 4)
  genSurvive <- lapply(genRepro, the.survivors, phase = pha, max.N = max.N)
  gen <- reduce(genSurvive, em.gl.join)
  gen@other$ind.metrics$generation <- i
  gen@other$ind.metrics$m <- m.ind
  gen@other$ind.metrics$Nmax <- max.N
  gen@other$ind.metrics$phaseNo <- icon$phaseNo
  gen@other$ind.metrics$phase <- icon$phase
  gen@other$ind.metrics$rep <- icon$rep
  
  sim[[i]] <- gen
  cat(paste("  Generation:", i, "\n"))

}
})

sapply(sim, nInd)/6


sapply(sim, nInd)*0.1
sim2 <- lapply(sim, function(x)  x[sample(1:nInd(x), ceiling(nInd(x)*0.2))])
system.time(fst <- lapply(sim2, gl.fst.pop))
system.time(fstall <- lapply(sim, gl.fst.pop))

sapply(fstall, function(x) mean(x$Fsts, na.rm = T))
#118 = 2mins for sim2

simgens <- reduce(sim, em.gl.join)
pop(simgens) <- simgens@other$ind.metrics$phaseNo

glPhase <- seppop(simgens)
