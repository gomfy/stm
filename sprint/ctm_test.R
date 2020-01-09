library(topicmodels)
library(stm)
library(stringr)

set.seed(2138)
heldout <- make.heldout(poliblog5k.docs, poliblog5k.voc)
slam <- convertCorpus(heldout$documents, heldout$vocab, type="slam")
control_CTM_VEM <- list(estimate.beta=TRUE, verbose=1, seed=as.integer(2138), nstart=1L, best=TRUE, var=list(iter.max=20, tol=1e-6), em=list(iter.max=1000, tol=1e-3), initialize="random", cg=list(iter.max = -1, tol = 1e-6))

k=143
local <- CTM(slam, k = k, control = control_CTM_VEM)
local_name <- str_c("ctm_blei_res_K", k, ".rda")
save(local, file = local_name)

# 
# cli = list()
# L <- 2
# U <- 10
# for (k in seq(L,U)) {
# 	cli <- CTM(slam, k = k, control = control_CTM_VEM)
# 	clin = str_c("ctm_blei_res_K", k, ".rda")
# 	save(cli, file = clin)
# }

