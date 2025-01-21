library(abind)
options(digits = 4)

#Set working directory
wd = ""
lambda = "1"
target = "mus"
setwd(wd)


#Load previously saved data and results
load(paste0(wd,"/snippets.RData"))
# load(paste0(wd,"/embeddings.RData"))
# load(paste0(wd,"/EDiSC_fit_data.RData"))

seed = ifelse(target == "bank", 100, 200)
load(paste0(wd,"/EDiSC_fit_seed_",seed,".RData"))


#############################################################################
# Extract Stan results ---------------------------------------------------- #
#############################################################################

# PERMUTED ----------------------------------------------------------------

#Extract posterior
EDiSC_posterior = rstan::extract(EDiSC_fit, permuted = TRUE)

#Rearrange posterior samples and variable names into familiar forms
num.samples = (EDiSC_fit@sim$iter - EDiSC_fit@sim$warmup) / EDiSC_fit@sim$thin
burn.in = 0

sense.probs.sim = array(data = aperm(exp(EDiSC_posterior$log_sense_probs), c(1,3,2)),
                        dim = c(num.samples, num.senses, num.snippets),
                        dimnames = list(Sample = 1:num.samples,
                                        Sense  = 1:num.senses, SnippetID = snippets$SnippetID))

psi.tilde.sim = array(data = aperm(EDiSC_posterior$psi_tilde, c(1,4,2,3)),
                      dim = c(num.samples, num.words, num.senses, num.periods),
                      dimnames = list(Sample = 1:num.samples,
                                      Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods))

phi.tilde.sim = array(data = aperm(EDiSC_posterior$phi_tilde, c(1,4,2,3)),
                      dim = c(num.samples, num.senses, num.genres, num.periods),
                      dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses,
                                      Genre = 1:num.genres, Time = 1:num.periods))

# chi.sim = array(data = aperm(EDiSC_posterior$chi, c(1,3,2)),
#                 dim = c(num.samples, embedding.dim, num.senses),
#                 dimnames = list(Sample = 1:num.samples, Embedding = 1:embedding.dim, Sense = 1:num.senses))
# 
# theta.sim = array(data = aperm(EDiSC_posterior$theta, c(1,3,2)),
#                   dim = c(num.samples, embedding.dim, num.periods),
#                   dimnames = list(Sample = 1:num.samples, Embedding = 1:embedding.dim, Time = 1:num.periods))
# 
# sigma.sim = array(data = EDiSC_posterior$sigma, dim = c(num.samples, num.words), 
#                   dimnames = list(Sample = 1:num.samples, Word = 1:num.words))
# 
# #Other diagnostics
# Run.Time = lubridate::as.duration(rstan::get_elapsed_time(EDiSC_fit)[2]) #run time for sampling phase only


# UNPERMUTED --------------------------------------------------------------

#Run above code then the following
sense.probs.sim[] = sense.probs.sim[order(EDiSC_fit@sim$permutation[[1]]),,]
psi.tilde.sim[] = psi.tilde.sim[order(EDiSC_fit@sim$permutation[[1]]),,,]
phi.tilde.sim[] = phi.tilde.sim[order(EDiSC_fit@sim$permutation[[1]]),,,]
# chi.sim[] = chi.sim[order(EDiSC_fit@sim$permutation[[1]]),,]
# theta.sim[] = theta.sim[order(EDiSC_fit@sim$permutation[[1]]),,]
# sigma.sim[] = sigma.sim[order(EDiSC_fit@sim$permutation[[1]]),]


# Clear memory ------------------------------------------------------------
rm(EDiSC_fit, EDiSC_posterior); gc()


#############################################################################
# Posterior means --------------------------------------------------------- #
#############################################################################

#posterior mean phi.tilde
phi.tilde.post.mean = sapply(1:num.periods, function(t) {
  sapply(1:num.genres, function(g) {apply(phi.tilde.sim[(burn.in+1):num.samples,,g,t], 2, mean)})}, simplify = "array" )
dimnames(phi.tilde.post.mean) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)

#posterior mean psi.tilde
psi.tilde.post.mean = sapply(1:num.periods, function(t) {
  sapply(1:num.senses, function(k) {apply(psi.tilde.sim[(burn.in+1):num.samples,,k,t], 2, mean)})}, simplify = "array" )
dimnames(psi.tilde.post.mean) = list(Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods)


#############################################################################
# Functions for PPCs ------------------------------------------------------ #
#############################################################################

#Memoise vectors of snippet time periods and genres
snippet.times  = snippets$Time
snippet.genres = snippets$genre


#Function to simulate z for one snippet given sense probabilities for that snippet
simulate.zd = function(sense.probs.d) {
  sample(1:num.senses, size = 1, prob = sense.probs.d)
}

#Function to simulate one z vector given sense probabilities for all snippets
simulate.z = function(sense.probs) {
  apply(sense.probs, 2, simulate.zd)
}


#function to simulate replicated data W given z and psi.tilde
simulate.data = function(z, psi.tilde) {
  
  #Sample words given z
  snippets.rep = lapply(1:num.snippets, function(d){sample(1:num.words, size = snippet.lengths[d], replace = TRUE, 
                                                           prob = psi.tilde[, z[d], snippet.times[d]])})
  snippets.rep = t(sapply(snippets.rep, "[", 1:snippet.length))
  
  #Return simulated snippets
  snippets.rep
  
}#simulate.data


# #Function to calculate log likelihood p(W,z|phi,psi) = p(W|z,psi)p(z|phi)
# get.log.lik = function(snippets, z, phi.tilde, psi.tilde) {
#   sum(sapply(1:num.snippets, function(d) {
#     sum(log(c(psi.tilde[as.numeric(snippets[d, 1:snippet.length]), z[d], snippet.times[d]], 
#               phi.tilde[z[d], snippet.genres[d], snippet.times[d]])), na.rm = TRUE)
#   }))
# }#get.log.lik


#Function to calculate log likelihood marginally over z: p(W|phi,psi) = Sum_z p(W|z,psi)p(z|phi)
get.log.lik.marg.over.z = function(snippets, phi.tilde, psi.tilde) {
  sum(sapply(1:num.snippets, function(d) {
    log(sum(apply(psi.tilde[as.numeric(snippets[d, 1:snippet.length]),,snippet.times[d]], 2, prod, na.rm = TRUE) *
              phi.tilde[,snippet.genres[d],snippet.times[d]]))
  }))
}#get.log.lik.marg.over.z


#############################################################################
# Simulate replicated data ------------------------------------------------ #
#############################################################################

#Simulate one z vector for each MCMC sample
z.sim = array(dim = c(num.samples, num.snippets), 
              dimnames = list(Sample = 1:num.samples, SnippetID = 1:num.snippets))
set.seed(100)
z.sim[] = t(apply(sense.probs.sim, 1, simulate.z))


#Simulate one set of replicated snippets for each MCMC sample
set.seed(200)
snippets.sim = aperm(sapply(1:num.samples, function(i) {simulate.data(z.sim[i,], psi.tilde.sim[i,,,])}, simplify = "array"), c(3,1,2))


#Save simulated data
save(z.sim, snippets.sim, file = paste0(wd,"/simulated.data.for.ppc.v4.RData"))
# load(paste0(wd,"/simulated.data.for.ppc.v4.RData"))


#############################################################################
# Evaluate PPC diagnostics ------------------------------------------------ #
#############################################################################

#Sub-sample (optional - quicker for testing)
# set.seed(500)
# samples = sample(1:num.samples,50)
samples = (burn.in+1):num.samples


#log likelihood at observed W, evaluated using phi.tilde and psi.tilde MCMC samples
log.lik.sim = colSums(log(apply(sense.probs.sim[samples,,], 1, colSums)))

#log likelihood at observed W, evaluated using posterior mean phi.tilde and psi.tilde
log.lik.mean = get.log.lik.marg.over.z(snippets, phi.tilde.post.mean, psi.tilde.post.mean)

#posterior predictive log likelihood at each replicated W, predicted and evaluated using using phi.tilde and psi.tilde MCMC samples
log.lik.ppc.pred.mcmc.eval.mcmc = sapply(samples, function(i) {
  get.log.lik.marg.over.z(snippets.sim[i,,], adrop(phi.tilde.sim[i,,,,drop=FALSE], drop=1), psi.tilde.sim[i,,,])
})


#Save log likelihoods 
save(log.lik.sim, log.lik.mean, log.lik.ppc.pred.mcmc.eval.mcmc, file = paste0(wd,"/log.lik.ppc.v4.RData"))
# load(paste0(wd,"/log.lik.ppc.v4.RData"))


#############################################################################
# Plots and p-values ------------------------------------------------------ #
#############################################################################

par(cex = 0.8)

pvalue = mean(log.lik.ppc.pred.mcmc.eval.mcmc < log.lik.mean)

plot.title = paste0(target, ", λ = ", lambda, ", p-value = ", round(pvalue,3))

dens.log.lik.sim = density(log.lik.sim)
dens.log.lik.ppc.pred.mcmc.eval.mcmc = density(log.lik.ppc.pred.mcmc.eval.mcmc)

plot(NULL, main = plot.title, xlab = "log likelihood", ylab = "density",
     xlim = range(c(log.lik.mean,
                    dens.log.lik.sim$x,
                    dens.log.lik.ppc.pred.mcmc.eval.mcmc$x
     )),
     ylim = range(c(dens.log.lik.sim$y, 
                    dens.log.lik.ppc.pred.mcmc.eval.mcmc$y
     )))

lines(dens.log.lik.sim, col = "black", lty = 1)
lines(dens.log.lik.ppc.pred.mcmc.eval.mcmc, col = "red", lty = 1)
abline(v = log.lik.mean, col = "black", lty = 2)

legend("topright",
  legend = c(
    "obs, at mcmc samples",
    "rep, at mcmc samples",
    "obs, at post mean"
  ),
  lty = c(1, 1, 2),
  col = c("black", "red", "black")
)

