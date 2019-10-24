#!/usr/bin/env Rscript

library(stm)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("Must pass the K parameter", call.=FALSE)
} 

K = as.numeric(args[1])
NUM_SAMPLES = 25 

print(sprintf('Deterministic: K = %d', K))
deter <- stm(poliblog5k.docs, poliblog5k.voc, K=K, #prevalence=~rating, 
             data=poliblog5k.meta, max.em.its = 100)

rand = list()
for(i in 1:NUM_SAMPLES)
{
    seed = sample(1:1e7,1)
    print(sprintf('Random: K = %d, sample %d of %d', K, i, NUM_SAMPLES))
    rand[[i]] <- stm(poliblog5k.docs, poliblog5k.voc, K=K, #prevalence=~rating, 
                     data=poliblog5k.meta, seed=seed, #seed=13, 
                     max.em.its = 100,
                     control=list(randomize=TRUE))
}


results = list(deter=deter, rand=rand)

saveRDS(results, sprintf('K_sweep_neumaier/K_sweep_%d.rds', K))
