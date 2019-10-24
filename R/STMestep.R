#
#E-Step for a Document Block
#[a relatively straightforward rewrite of previous
# code with a focus on avoiding unnecessary computation.]

#Input: Documents and Key Global Parameters
#Output: Sufficient Statistics

# Approach:
# First we pre-allocate memory, and precalculate where possible.
# Then for each document we:
#  (1) get document-specific priors, 
#  (2) infer doc parameters, 
#  (3) update global sufficient statistics
# Then the sufficient statistics are returned.

#Let's start by assuming its one beta and we may have arbitrarily subset the number of docs.
#' @importFrom PreciseSums fsum
estep <- function(documents, beta.index, update.mu, #null allows for intercept only model  
                       beta, lambda.old, mu, sigma, 
                       order_sigma, order_beta, randomize,
                       verbose) {
  #quickly define useful constants
  V <- ncol(beta[[1]])
  K <- nrow(beta[[1]])
  N <- length(documents)
  A <- length(beta)
  ctevery <- ifelse(N>100, floor(N/100), 1)
  if(!update.mu) mu.i <- as.numeric(mu)
  # 1) Initialize Sufficient Statistics 
  if(order_sigma) {
    sigma.ss <- n_mat_sum(diag(0, nrow=(K-1)))
  } else {
    sigma.ss <- diag(0, nrow=(K-1))
  }
  if(order_beta) {
    beta.ss <- vector(mode="list", length=A)
    for(i in 1:A) {
      beta.ss[[i]] <- n_mat_sum(matrix(0, nrow=K,ncol=V))
    }
  } else {
    beta.ss <- vector(mode="list", length=A)
    for(i in 1:A) {
      beta.ss[[i]] <- matrix(0, nrow=K,ncol=V)
    }
  }
  bound <- vector(length=N)
  lambda <- vector("list", length=N)
  
  # 2) Precalculate common components
  sigobj <- try(chol.default(sigma), silent=TRUE)
  if(class(sigobj)=="try-error") {
    sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    siginv <- solve(sigma)
  } else {
    sigmaentropy <- sum(log(diag(sigobj)))
    siginv <- chol2inv(sigobj)
  }
  # 3) Document Scheduling
  # For right now we are just doing everything in serial.
  # the challenge with multicore is efficient scheduling while
  # maintaining a small dimension for the sufficient statistics.
  if(randomize) {
    vec <- sample(1:N, N) 
  } else {
    vec <- 1:N
  }
  
  # We need to keep track of all components of the summations. Probably can find 
  # an online algorithm somewhere but this is what for now.
  betas = array(rep(0, N*prod(dim(beta.ss[[1]]))), c(dim(beta.ss[[1]]), N))
  sigmas = array(rep(0, N*prod(dim(sigma.ss))), c(dim(sigma.ss), N))
  
  for(i in vec) {
    #update components
    doc <- documents[[i]]
    words <- doc[1,]
    aspect <- beta.index[i]
    init <- lambda.old[i,]
    if(update.mu) mu.i <- mu[,i]
    beta.i <- beta[[aspect]][,words,drop=FALSE]
    
    #infer the document
    doc.results <- logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, 
                                  doc=doc, sigmaentropy=sigmaentropy)
    
    # update sufficient statistics 
    if(order_sigma) {
      sigma.ss <- n_mat_sum(sigma.ss[[1]], sigma.ss[[2]], doc.results$eta$nu)
    } else {
      #sigma.ss <- sigma.ss + doc.results$eta$nu
      sigmas[1:dim(sigma.ss)[1], 1:dim(sigma.ss)[2], i] <- doc.results$eta$nu
    }
    if(order_beta) {
      #more efficient than this would be to stack all the C's underneath
      #betas
      o_beta <- n_mat_sum(beta.ss[[aspect]][[1]][,words], 
                          beta.ss[[aspect]][[2]][,words], 
                          doc.results$phis)
      beta.ss[[aspect]][[1]][,words] <- o_beta[[1]]
      beta.ss[[aspect]][[2]][,words] <- o_beta[[2]]
    } else {
      #beta.ss[[aspect]][,words] <- doc.results$phis + beta.ss[[aspect]][,words]
      betas[1:dim(betas)[1], words, i] <- doc.results$phis
    }
    bound[i] <- doc.results$bound
    lambda[[i]] <- c(doc.results$eta$lambda)
    if(verbose && i%%ctevery==0) cat(".")
  }
  if(verbose) cat("\n") #add a line break for the next message.
  
  sigma_exact = matrix(0, dim(sigma.ss)[1], dim(sigma.ss)[2])
  for(i in 1:dim(sigma.ss)[1])
    for(j in 1:dim(sigma.ss)[2])
      sigma_exact[i,j] = PreciseSums::neumaierSum(sigmas[i, j, vec])
  
  beta_exact = matrix(0, dim(beta.ss[[aspect]])[1], dim(beta.ss[[aspect]])[2])
  for(i in 1:dim(beta_exact)[1])
    for(j in 1:dim(beta_exact)[2])
      beta_exact[i,j] = PreciseSums::neumaierSum(betas[i, j, vec])
  
  sigma.ss <- sigma_exact
  beta.ss[[aspect]] <- beta_exact
  
  #4) Combine and Return Sufficient Statistics
  lambda <- do.call(rbind, lambda)
  return(list(sigma=sigma.ss, beta=beta.ss, bound=bound, lambda=lambda,
              vec=vec))
}
