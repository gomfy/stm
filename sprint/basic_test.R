library(stm)
detero <- stm(poliblog5k.docs, poliblog5k.voc, K=5, prevalence=~rating, data=poliblog5k.meta, max.em.its = 100, control=list(method="BFGS", order_sigma=TRUE, order_beta=TRUE))
rando <- stm(poliblog5k.docs, poliblog5k.voc, K=5, prevalence=~rating, data=poliblog5k.meta, seed=13, max.em.its = 100, control=list(method="BFGS", randomize=TRUE, order_sigma=TRUE, order_beta=TRUE))
all.equal(detero, rando)

deter <- stm(poliblog5k.docs, poliblog5k.voc, K=5, prevalence=~rating, data=poliblog5k.meta, max.em.its = 100, control=list(method="BFGS"))
rand <- stm(poliblog5k.docs, poliblog5k.voc, K=5, prevalence=~rating, data=poliblog5k.meta, seed=13, max.em.its = 100, control=list(method="BFGS", randomize=TRUE))
all.equal(deter, rand)