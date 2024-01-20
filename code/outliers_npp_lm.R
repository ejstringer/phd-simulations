# grids sampled ------------ 
source('./code/libraries.R')
source('../phd-analysis2/code/functions_genetic-data-creation.R')
library(fuzzySim)
ph <- em.filtering('sy') 
ph@other$ind.metrics$trip <- ymd(ph@other$ind.metrics$trip)
ph@other$ind.metrics<- ph@other$ind.metrics %>% left_join(em.rain_caps_phase())

ph@other$ind.metrics$phaseNo <- factor(ph@other$ind.metrics$phaseNo,
                                       levels = paste0(rep(c('L', 'I', 'D'), 3),
                                                       rep(1:3, each = 3)))

npp <- read.csv('../phd-analysis2/output/npp_means.csv')
npp

ph@other$ind.metrics <- ph@other$ind.metrics %>% left_join(npp)
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

# filter real data ----------------------
ph 
gridNames %>% length
glBase@loc.names

ph <- gl.keep.loc(ph, loc.list = glBase@loc.names)
ph_real <- gl.keep.pop(ph, pop.list = gridNames, as.pop = 'gridId')
nInd(ph)-nInd(ph_real)

# sep phase --------------

pop(ph_real) <- ph_real@other$ind.metrics$phaseNo
ph_phase <- seppop(ph_real) 
ph_phase <- lapply(ph_phase, em.gl_assign_pop)

# phase summary
meta <- ph@other$ind.metrics %>% 
  group_by(period, phaseNo, phase, npp) %>% 
  summarise(n = n(),
            ngrids = length(unique(gridId)))

# all freq ----------

alfphase <- pblapply(ph_phase, gl.alf)

for(i in 1:9){
  
  a <- alfphase[[i]]
  alfphase[[i]]$loci <- rownames(a)
  alfphase[[i]]$phaseNo <- names(alfphase)[i]
  
}

alfdf <- do.call('bind_rows', alfphase) %>% 
  left_join(meta) %>% 
  mutate(npp.log = log(npp))
rownames(alfdf) <- NULL
table(alfdf$phaseNo)
alfdf$phaseNo <- factor(alfdf$phaseNo,
                        levels = levels(ph@other$ind.metrics$phaseNo))



# models ----------
em.alf_npp_model <- function(a1){
  #a1 <- alfdf %>% filter(loci == alfdf$loci[2])
  
  m <- lm(alf2 ~ npp.log, a1)
  m_summary <- summary(m)
  
  data.frame(loci = a1$loci[1], R2 = m_summary$r.squared,
             intercept = coef(m_summary)[1,1], slope = coef(m_summary)[2,1], 
             df = m_summary$df[2], tvalue =  coef(m_summary)[2,3],
             pvalue =  coef(m_summary)[2,4])
  
}


alfloci <- split(alfdf, alfdf$loci)

library(fuzzySim)
my_pvalues <- data.frame(var = letters[1:5], pval = c(0.02, 0.004, 0.07, 0.03, 0.05))
FDR(pvalues = my_pvalues)

m.alfnpp <- lapply(alfloci, em.alf_npp_model) %>% 
  do.call('rbind', .)
1210/24388

m.alfnpp %>% filter(pvalue < 0.05) %>% nrow
m.alfnpp$pvalue %>% hist
m.alfnpp %>% filter(intercept > 0.1 , intercept < 0.9, pvalue < 0.001)

test_pvalues <- m.alfnpp[,c('loci', 'pvalue')]
padjusted <- FDR(pvalues = test_pvalues)
padjusted$exclude$p.adjusted %>% 
  table %>% data.frame(p = as.numeric(names(.), f = .)) %>% 
  ggplot(aes(p, Freq)) +
  geom_point()

notOutliers <- padjusted$exclude[padjusted$exclude$p.adjusted < 0.7,]
# plot --------
phaseCol <- c(low = terrain.colors(10)[7],
              increase = terrain.colors(10)[1],
              decrease = terrain.colors(10)[4])


em.alf_npp_model(filter(alfdf, loci == '15083484-51-G/A', phaseNo != 'L1'))
a1 <- filter(alfdf, loci %in% c('50657404-35-C/T','15083484-51-G/A'))
ggplot(a1, aes(npp, alf2, colour = loci, group = loci)) +
  geom_smooth(method = 'lm', colour = 'grey50', se = T) +
  geom_hline(yintercept = 0, colour = 'grey', lty = 2)+
  geom_point(size = 2.5)+
  scale_x_log10()+
  ylim(0, 0.08)+
  #scale_colour_manual(values = phaseCol)+
  theme_classic()+
  theme(legend.position = 'bottom')+
  
  ggplot(a1, aes(phaseNo, alf2, colour = phase, group = loci)) +
  #geom_smooth(method = 'lm', colour = 'grey50', se = T) +
  geom_line()+
  ylim(0, 0.08)+
  geom_hline(yintercept = 0, colour = 'grey', lty = 2)+
  geom_point(aes(size = n))+
  scale_colour_manual(values = phaseCol)+
  theme_classic()+
  theme(axis.title.y = element_blank())

ph@other$loc.metrics %>% filter(uid %in%  c('15083484-51', '50657404-35'))

lookout <- sub('-C/T', '', rownames(notOutliers))#[-c(4,8,10)])
lookout <- sub('-G/C', '', lookout)
lookout <- sub('-G/A', '', lookout)
lookout <- sub('-C/A', '', lookout)
ph@other$loc.metrics %>%
  filter(uid %in% lookout) %>%
  .[,c(1,9:10)]
ggsave('./figures/outliers.png', units = 'cm', height = 12, width = 18)

## check ----
a1 <- filter(alfdf, loci %in% c('50657404-35-C/T','15083484-51-G/A', '15083908-17-G/A'))
a1 <- filter(alfdf, loci %in% rownames(notOutliers)[-c(4,8,10)])
a1 <- filter(alfdf, loci %in% rownames(notOutliers))
ggplot(a1, aes(npp, alf2, colour = loci, group = loci)) +
  geom_smooth(method = 'lm', colour = 'grey50', se = T) +
  #geom_hline(yintercept = 0, colour = 'grey', lty = 2)+
  geom_point(size = 3)+
  scale_x_log10()+
  theme_classic()+
  facet_wrap(~loci, scale = 'free')+
  theme(legend.position = 'bottom')


ggplot(a1, aes(phaseNo, alf2, colour = phase, group = loci)) +
  #geom_smooth(method = 'lm', colour = 'grey50', se = T) +
  geom_line(size = 1)+
  geom_hline(yintercept = 0, colour = 'grey', lty = 2)+
  geom_point(aes(size = n))+
  scale_colour_manual(values = phaseCol)+
  theme_classic()+
  facet_wrap(~loci, scale = 'free')+
  theme(axis.title.y = element_blank())

ggplot(m.alfnpp, aes(tvalue, R2))+
  geom_point()+
  theme_classic()


# slopes ---------

m.alfnpp$slope %>% abs %>% hist
m.alfnpp %>% filter(slope > 0.04, pvalue < 0.05) %>% 
  arrange(slope)

a1 <- filter(alfdf, loci %in% c('15093965-10-A/G', '15089213-50-A/G'))
ggplot(a1, aes(npp, alf2, colour = loci, group = loci)) +
  geom_smooth(method = 'lm', colour = 'grey50', se = T) +
  geom_hline(yintercept = 0, colour = 'grey', lty = 2)+
  geom_point(size = 2.5)+
  scale_x_log10()+
  #scale_colour_manual(values = phaseCol)+
  theme_classic()+
  theme(legend.position = 'bottom')+
  
  ggplot(a1, aes(phaseNo, alf2, colour = phase, group = loci)) +
  #geom_smooth(method = 'lm', colour = 'grey50', se = T) +
  geom_line()+
  geom_hline(yintercept = 0, colour = 'grey', lty = 2)+
  geom_point(aes(size = n))+
  scale_colour_manual(values = phaseCol)+
  theme_classic()+
  theme(axis.title.y = element_blank())


lookout <- sub('-A/G', '', a1$loci)#[-c(4,8,10)])
ph@other$loc.metrics %>%
  filter(uid %in% lookout)

# r sqaured ----
m.alfnpp$R2 %>% hist
topR2 <- m.alfnpp %>% filter(R2 > 0.8) %>% 
  arrange(desc(R2))
topR2 %>% nrow
topR2 %>% 
  ggplot(aes(pvalue, R2)) +
  geom_point()+
  theme_classic()

a1 <- filter(alfdf, loci %in% topR2$loci)
a1 <- filter(alfdf, loci %in% outsScaf)

ggplot(a1, aes(npp, alf2, colour = loci, group = loci)) +
  geom_smooth(method = 'lm', colour = 'grey50', se = T) +
  geom_hline(yintercept = 0, colour = 'grey', lty = 2)+
  geom_point(size = 2.5)+
  scale_x_log10()+
  #scale_colour_manual(values = phaseCol)+
  theme_classic()+
  theme(legend.position = 'bottom')+
  
  ggplot(a1, aes(phaseNo, alf2, colour = phase, group = loci)) +
  #geom_smooth(method = 'lm', colour = 'grey50', se = T) +
  geom_line()+
  geom_hline(yintercept = 0, colour = 'grey', lty = 2)+
  geom_point(aes(size = n))+
  scale_colour_manual(values = phaseCol)+
  theme_classic()+
  theme(axis.title.y = element_blank())

lookout <- sub('-C/T', '',topR2$loci)#[-c(4,8,10)])
lookout <- sub('-G/C', '', lookout)
lookout <- sub('-G/A', '', lookout)
lookout <- sub('-C/A', '', lookout)
lookout <- sub('-A/G', '', lookout)
lookout <- sub('-G/T', '', lookout)
ph@other$loc.metrics %>%
  filter(uid %in% lookout) %>% 
  arrange(Chrom_Pseudomys_desertor_wtdbg2)-> outs
  outs$Chrom_Pseudomys_desertor_wtdbg2 %>% factor %>%  table %>% 
  names
outs 
outs$uid[c(11:14,19)]
outs$uid[16:18]
  ph@other$loc.metrics$Chrom_Pseudomys_desertor_wtdbg2 %>% table

  outsScaf <- topR2$loci[sapply(outs$uid[c(11,19)], function(x) grep(x, topR2$loci))]
  outsScaf <- topR2$loci[sapply(outs$uid[12:14], function(x) grep(x, topR2$loci))]
 outsScaf <- topR2$loci[sapply(outs$uid[16:18], function(x) grep(x, topR2$loci))]
  
 
 
 # syoung --------
 m.alfnpp$R2 %>% hist
 topR2 <- m.alfnpp %>% filter(R2 > 0.7) %>% 
   arrange(desc(R2))
 topR2 %>% nrow
 topR2 %>% 
   ggplot(aes(pvalue, R2)) +
   geom_point()+
   theme_classic()
 
 a1 <- filter(alfdf, loci %in% topR2$loci)
 ggplot(a1, aes(npp, alf2, colour = loci, group = loci)) +
   geom_smooth(method = 'lm', colour = 'grey50', se = T) +
   geom_hline(yintercept = 0, colour = 'grey', lty = 2)+
   geom_point(size = 2.5)+
   scale_x_log10()+
   #scale_colour_manual(values = phaseCol)+
   theme_classic()+
   theme(legend.position = 'bottom')+
   
   ggplot(a1, aes(phaseNo, alf2, colour = phase, group = loci)) +
   #geom_smooth(method = 'lm', colour = 'grey50', se = T) +
   geom_line()+
   geom_hline(yintercept = 0, colour = 'grey', lty = 2)+
   geom_point(aes(size = n))+
   scale_colour_manual(values = phaseCol)+
   theme_classic()+
   theme(axis.title.y = element_blank())
 
 lookout <- sub('-C/G', '',topR2$loci)#[-c(4,8,10)])
 lookout <- sub('-T/A', '', lookout)
 lookout <- sub('-G/A', '', lookout)

 phsy@other$loc.metrics[phsy@other$loc.metrics$uid %in% lookout,c(1,24,5,6)]
 
 ggplot(a1, aes(npp, alf2, colour = loci, group = loci)) +
   geom_smooth(method = 'lm', colour = 'grey50', se = T) +
   #geom_hline(yintercept = 0, colour = 'grey', lty = 2)+
   geom_point(size = 2.5)+
   scale_x_log10()+
   facet_wrap(~loci, scale = 'free')+
   #scale_colour_manual(values = phaseCol)+
   theme_classic()+
   theme(legend.position = 'bottom')

 m.alfnpp$slope %>% hist 
 