#load file
getwd()
source('./code/libraries_dungog.R')
source('./code/functions_simulate_geneflow.R')
ph <- readRDS('./output/syoung_filtered_genotypes_phases.rds')

ph@other$ind.metrics <- ph@other$ind.metrics %>% 
  mutate(phaseNo = factor(phaseNo,
                          levels  = factor(paste0(rep(c('L', 'I', 'D'), 3),
                                                  rep(1:3, each = 3)))))

em.mRate <- function(fst, ne) (1/fst-1)/(4*ne)

em.gl_assign_pop <- function (glx, define.pop = "gridId") 
{
  if (length(define.pop) == 1) {
    if (!define.pop %in% names(glx@other$ind.metrics)) {
      stop(paste(define.pop, "not in ind.metrics; names(glx@other$ind.metrics)"))
    }
    pop(glx) <- glx@other$ind.metrics[, define.pop]
  }
  else {
    pop(glx) <- define.pop
  }
  return(glx)
}

