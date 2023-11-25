# functions -------
em.simulation.gdist <- function(zero, mat, N_values, m_values){
  sim <- list()
  x <- 1
  m.iteration <- 0
  system.time({
    for(N in N_values){
      cat(paste("Population size:", N, "\n"))
      for(m in m_values){
        m.ind <- m
        gen <- zero #gene.flow(zero, mat, m.ind)
        m.iteration <- m.iteration + 1
        cat(paste(" Migraton:", m, "\n"))
        for(i in 1:5){
          gen <- gene.flow(gen, mat, m.ind)
          pops <- seppop(gen)
          genRepro <- lapply(pops, next.gen, N.size = N)
          gen <- reduce(genRepro, em.gl.join)
          gen@other$ind.metrics$generation <- i
          gen@other$ind.metrics$m <- m.ind
          gen@other$ind.metrics$N <- N
          gen@other$ind.metrics$iteration <- m.iteration
          sim[[x]] <- gen
          if(i == 5) cat(paste("  Generation:", i, "\n"))
          x <- x+1
        }
      }
    }
  })
  cat(length(sim))
  return(sim)
}

em.gdist <- function(simx){
  df <- matrixConvert(dist(as.matrix(simx)-1),
                      colname = c('id1','id2','gdist'))
  df$m <- simx@other$ind.metrics$m[1]
  df$gen <- simx@other$ind.metrics$gen[1]
  df$N <- simx@other$ind.metrics$N[1]
  df$iteration <- simx@other$ind.metrics$iteration[1]
  return(df)
}
# library -------------
source('./code/libraries.R')
source('./code/functions_genetic-data-creation.R')
source('./code/functions-npp-data-creation.R')
source('./code/functions_simulate_geneflow.R')

# genetic --------------
ph <- em.filtering('ph')

ph2 <- gl.keep.pop(ph, as.pop = "trip", pop.list = "2007-09-01")
pop(ph2) <- ph2@other$ind.metrics$gridId
ph2 <- gl.filter.callrate(ph2, threshold = 1, plot = F)
glBase <- gl.keep.pop(ph2, as.pop = "sex", pop.list = c("m", "f"))
table(glBase@other$ind.metrics$sex, exclude = F)

# mat ---------------------
grids <- read.csv("../data/em_gridcoortable.csv")
grids[grids$coor== "site mean",c("lat", "lon")] <- NA
xyGrids <- rgdal::project(as.matrix(grids[,c("lon", "lat")]),
                          proj = '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
rownames(xyGrids) <- grids$gridId
#xy <- xyGrids[rownames(xyGrids) %in% c("KSE2", "MC4",  "SS3",  "FRN1", "FRS2", "MC2"),]
xy <- xyGrids[rownames(xyGrids) %in% c("MC2", "MC4"),]
nrow(xy)

lcost <- xy
mat <- as.matrix(dist(lcost))
mat[mat != 0] <- 1/(nrow(mat)-1)

colSums(mat)
mat
## zero ----------

sites <- rownames(xy) # populations
p <- lapply(sites, function(x) gen.zero(glBase, x))
zero <- reduce(p, em.gl.join)
table(pop(zero)) 

# Scenario 1 ----------
sim <- em.simulation.gdist(zero,mat, 
                    N_values = rep(c(10,100), each = 4),
                    m_values = seq(0,0.5, 0.1))

sim <- em.simulation.gdist(zero,mat, 
                           N_values = rep(c(10,100), each = 4),
                           m_values = seq(0,0.5, 0.1))

system.time(gdist <- lapply(sim, em.gdist))

gdist <- readRDS('./output/gdist_simulation_geneflow_1st_N10_50_100_2sites.rds')
  dfgdist <- gdist %>% 
  do.call('rbind', .) %>%
  mutate(site1 = str_sub(id1, 1,3),
         site2 = str_sub(id2, 1,3),
         comparison = factor(ifelse(site1!=site2, 'between', 'within'),
                             levels= c('within', 'between')))

names(dfgdist)
mean.gdist <- dfgdist$gdist %>% mean()
### gdist x generation --------

df2 <- dfgdist %>%  
  group_by(gen, m, N,comparison) %>% 
  summarise(dist = mean(gdist),
            std = sd(gdist),
            n = n(),
            se = std/sqrt(n),
            lower = dist-se*1.96,
            upper = dist+se*1.96) 
ggplot(df2, aes(gen, dist, fill = comparison))+
  #geom_errorbar(aes(ymin=lower, ymax = upper), width = 0)+
  geom_line()+
  geom_point(pch = 21, size = 3)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c('grey90', '#F8766D'))+
  facet_grid(N~m, scale = 'free_y')+
  scale_y_log10()

### gdist x iteration ---------

df3 <- dfgdist %>% 
  group_by(iteration, m, N,comparison) %>% 
  summarise(dist = mean(gdist),
            skew = moments::skewness(gdist))
ggplot(df3, 
       aes(m, dist, fill = factor(m),
           colour = factor(N)), group = m)+
  geom_smooth(method = 'loess', aes(group = comparison), 
              fill = 'grey', se = T, colour = 'grey50')+
  geom_boxplot(colour = 'black', width = 0.05)+
  # #geom_line(size = 1)+
  geom_hline(yintercept = mean.gdist,
             linetype = 'dashed', colour ='#F8766D')+
  #geom_point(pch = 21, size = 3, colour = 'black')+
  theme_bw()+
  #scale_x_log10()+
  facet_grid(N~comparison, scale = 'free')+
  theme(panel.grid = element_blank())

ggplot(filter(df3, m > 0), 
       aes(m, skew, fill = factor(m),
           colour = factor(m), group = m))+
  geom_smooth(method = 'lm', aes(group = comparison), 
              fill = 'grey', se = T, colour = 'grey50')+
  geom_point(pch = 21, size = 2)+
  theme_bw()+
  geom_boxplot(colour = 'black', width = 0.025,position = position_dodge(0.75))+
  # scale_colour_manual(values = c('orange', 'forestgreen'))+
  # scale_fill_manual(values = c('honeydew', 'lightgreen'))+
  #scale_x_log10()+
  facet_grid(N~comparison, scale = 'free_y')+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

### gdist x pop size ----------
df4 <- dfgdist %>% 
  group_by(m, N,comparison) %>% 
  summarise(dist = mean(gdist),
            skew = moments::skewness(gdist))

ggplot(df4, 
       aes(N, dist, fill = factor(m),
           colour = factor(m), group = m))+
  geom_smooth(method = 'loess', se = F)+
  geom_hline(yintercept = mean.gdist,linetype = 'dashed', colour ='#F8766D')+
  geom_point(pch = 21, size = 3, colour = 'black')+
  theme_bw()+
  scale_x_log10()+
  facet_grid(~comparison)+
  theme(panel.grid = element_blank())


ggplot(filter(df4, m < 0.3 & m > 0, N != 50), 
       aes(N, dist, fill = factor(m),
           colour = factor(m), group = m))+
  geom_smooth(method = 'lm', se = F)+
  geom_hline(yintercept = mean.gdist,linetype = 'dashed', colour ='#F8766D')+
  geom_point(pch = 21, size = 3, colour = 'black')+
  theme_bw()+
  scale_x_log10()+
  facet_grid(~comparison)+
  theme(panel.grid = element_blank())

ggplot(filter(df4, m > -1), 
       aes(m, dist, fill = comparison,
           colour = comparison, group = comparison))+
  geom_smooth(method = 'loess', se = T)+
  geom_point(pch = 21, size = 3, colour = 'black')+
  theme_bw()+
  scale_colour_manual(values = c('grey70', '#F8766D'))+
  scale_fill_manual(values = c('grey70', '#F8766D'))+
  facet_wrap(~N, scale = 'free_y')+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())


ggplot(df4, 
       aes(N, skew, fill = factor(m),
           colour = factor(m), group = m))+
  geom_smooth(method = 'lm', se = F)+
  geom_hline(yintercept = 0,linetype = 'dashed', colour ='grey')+
  geom_point(pch = 21, size = 3, colour = 'black')+
  theme_bw()+
  scale_x_log10()+
  facet_grid(~comparison)+
  theme(panel.grid = element_blank())

hex <- scales::hue_pal()(5)

ggplot(filter(df4, m %in% c(0, 0.1, 0.5)), 
       aes(N, skew, fill = factor(m),
           colour = factor(m), group = m))+
  geom_smooth(method = 'lm', se = F)+
  geom_hline(yintercept = 0,linetype = 'dashed', colour ='grey')+
  geom_point(pch = 21, size = 3, colour = 'black')+
  theme_bw()+
  scale_x_log10()+
  scale_colour_manual(values = hex[c(1,2,5)])+
  scale_fill_manual(values = hex[c(1,2,5)])+
  facet_grid(~comparison)+
  theme(panel.grid = element_blank())
  