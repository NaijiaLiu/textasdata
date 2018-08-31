#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


#' @title Calculating Marginal Likelihiood or Topic Model 
#' @description Calculate the marginal likelihood.
#' @name Marginalikelihood
#' @author Naijia Liu
#' @usage null
#' @param wordspertopic  Number of words per topic (integer)
#' @param LDAresult      Results of LDA model
#' @param dtmatrix       Document term matrix
#' @param draws          Number of draws in MCMC
#' @param seed           Seed number
#' @param verbose        Print out results for each document if TRUE
#' @import topicmodels
#' @import MCMCpack
#' 
#'
#' @format null
#' @return
#' Return a matrix of marginal probability
#' Return a matrix of marginal probability (with penalization of unlikely documents)
#' Return a visualization of the best K
#' @examples
#' library(topicmodels)
#' library(MCMCpack)
#' data("AssociatedPress")
#' model <- LDA(AssociatedPress[1:10,],k=5, control = list(alpha=0.1, seed=99) ,  method = 'VEM')
#' Marginalikelihood(10, model, AssociatedPress[1:10,],1000,99,verbose = FALSE)
#'
#'
#' @examples
#' with raw text data
#' text <- lda::lexicalize(text)
#' dtm <- topicmodels::ldaformat2dtm(text$document, text$vocab)
#' model <- LDA(dtm, k=5, control = list(alpha=0.1, seed=99) ,  method = 'VEM')
#' Marginalikelihood(10, model, dtm, 1000, 99,verbose = FALSE)
#' 
#' @examples
#' 
#' @export
#'
#'
#'
Marginalikelihood <- function(wordspertopic, LDAresult ,dtmatrix, draws, seed, verbose){
  topicprob <- NULL
  Marginal <- NULL
  summation <- NULL
  for(d in 1:(LDAresult@Dim[1])){
    alpha <- LDAresult@call$control$alpha
    method <- LDAresult@call$method
    term1 <- ifelse(method=='VEM',LDAresult@loglikelihood[d],
                    LDAresult@loglikelihood/LDAresult@Dim[1])#loglik at that point
    term2a <- gamma(alpha*LDAresult@k) / (gamma(alpha))^LDAresult@k
    theta <- 1/LDAresult@k
    term2b <- (theta^(alpha-1))^LDAresult@k
    term2c <- (theta*(1/dtmatrix$ncol))^wordspertopic
    term2c <-( LDAresult@k *term2c ) #^dtmatrix$ncol
    term2d <- LDAresult@k*(term2b *term2c)
    term2 <- term2a * term2d
    #term2 <- sum(dtmatrix[d,]$i)/dtmatrix$ncol ## wrong way?
    term2 <- log(term2)
    
    set.seed(seed)
    gamma <- MCmultinomdirichlet(LDAresult@gamma[d,], alpha0 = rep(LDAresult@alpha, LDAresult@k), mc=draws, seed=seed)
    #  gamma <- rdirichlet(n = draws, alpha = LDAresult@gamma[d,])
    identify <- seq(1, draws, 5)
    newtheta <- colSums(gamma[identify,])/length(identify)  #Topicweights
   
    ##############################
  
    term3a <- gamma((LDAresult@alpha)*(LDAresult@k)) /(gamma(LDAresult@alpha))^LDAresult@k
    term3b <- prod(newtheta^(LDAresult@alpha-1))
    thetabeta <- NULL
    term3c <- NULL
    for(topic in 1:LDAresult@k) {
      testvalue <- sort(LDAresult@beta[topic,], decreasing = TRUE)[wordspertopic]
      #summary(LDAresult@beta[topic,])[5]
      set.seed(seed)
      betadraws <- MCmultinomdirichlet(exp(LDAresult@beta[topic,]), alpha0 = rep(1, LDAresult@Dim[2]), mc=draws, seed=seed)
      newbetas <- colSums(betadraws[identify,])/length(identify)
      
      for(word in 1:dtmatrix$ncol){
        thetabeta[word] <- ifelse(log(newbetas[word]) >= testvalue, 
                                  #  LDAresult@gamma[d,topic]* exp(LDAresult@beta[topic, word]),
                                  newtheta[topic] * newbetas[word],
                                  1)
      }
      term3c[topic] <- prod(thetabeta)
    }
    term3c <- sum(term3c)
    term3 <- term3a*term3b*term3c
    term3 <- log(term3)
    #########################
    
    
    Marginal[d] <- ifelse(term3 > -Inf, term1+term2-term3, 0 )
    if(verbose=="TRUE") print(c(term1, term2, term3a, term3b, term3c, term3, Marginal[d]))
  }
  #MarginalP <- sum(Marginal) 
  MarginalP <- sum(Marginal)/LDAresult@Dim[1]
  return(MarginalP)
}

##########################################################################
#######################Visulization#######################################
#' 
#' @title Visualizing Marginal Likelihiood for Topic Model
#' @name visualizeTMM
#' @author Naijia Liu
#' @usage null
#' @description This function visualizes the marginal likelihood value, and returns a vector of log-likelihood.
#' @param wordspertopic  Number of words per topic (integer)
#' @param dtmatrix       Document term matrix
#' @param draws          Number of draws in MCMC
#' @param seed           Seed number
#' @param verbose        Print out results for each document if TRUE
#' @param Kstart        starting value of topic number
#' @param Kend          ending value of topic number
#' @import topicmodels
#' @import MCMCpack
#' @import ggplot2
#' 
#' 
#' @examples
#' library(topicmodels)
#' library(MCMCpack)
#' data("AssociatedPress")
#' visualizeTMM(2, 10, dtmatrix = AssociatedPress[1:5,], 99, 5, 1000, FALSE)
#'
#' 
#' @export
#' 
visualizeTMM <- function(Kstart, Kend, dtmatrix, seed, wordspertopic, draws, verbose){
  index <- 0
  drawed <- NULL
  for(kk in seq(Kstart,Kend, 1)){
    index <- index +1
    
    simulation <- LDA(dtmatrix, control = list(alpha=0.1, seed=seed) , 
                      k = kk, method = 'VEM')
    drawed[index] <- Marginalikelihood(wordspertopic = wordspertopic, 
                                       LDAresult =  simulation, 
                                       dtmatrix = dtmatrix, 
                                       draws = draws, seed = seed, verbose = FALSE)
    if(verbose=="TRUE") print(paste('done with', index, sep = ' '))
  }
  
  PLOTing <- plot(drawed,  xlab = paste('Topic Number from ',Kstart,'to', Kend, sep = ' '),
                  ylab='Marginal Likelihood')
  
  return(drawed)
  print(PLOTing)
}
## attr(visualizeTMM, "help") <- 











