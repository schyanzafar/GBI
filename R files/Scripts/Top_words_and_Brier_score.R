library(tidyr)
library(plyr)
library(measures)
options(digits = 4)


#Load previously saved data and results
load("snippets.RData")
# load("embeddings.RData")
# load("EDiSC_fit_data.RData")


seeds = c(100,200)

for (seed in seeds) {

load(paste0("EDiSC_fit_seed_",seed,".RData"))

print(getwd())
print(paste("seed", seed))


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

# phi.tilde.sim = array(data = aperm(EDiSC_posterior$phi_tilde, c(1,4,2,3)),
#                       dim = c(num.samples, num.senses, num.genres, num.periods),
#                       dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses,
#                                       Genre = 1:num.genres, Time = 1:num.periods))
# 
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
# phi.tilde.sim[] = phi.tilde.sim[order(EDiSC_fit@sim$permutation[[1]]),,,]
# chi.sim[] = chi.sim[order(EDiSC_fit@sim$permutation[[1]]),,]
# theta.sim[] = theta.sim[order(EDiSC_fit@sim$permutation[[1]]),,]
# sigma.sim[] = sigma.sim[order(EDiSC_fit@sim$permutation[[1]]),]

# Clear memory ------------------------------------------------------------
rm(EDiSC_fit, EDiSC_posterior); invisible(gc())


#############################################################################
# psi and phi posterior means --------------------------------------------- #
#############################################################################

psi.tilde.post.mean = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k) 
  {apply(psi.tilde.sim[(burn.in+1):num.samples,,k,t], 2, mean)})}, simplify = "array" )
dimnames(psi.tilde.post.mean) = list(Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods)

# phi.tilde.post.mean = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) 
#   {apply(phi.tilde.sim[(burn.in+1):num.samples,,g,t], 2, mean)})}, simplify = "array" )
# dimnames(phi.tilde.post.mean) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)

# rm(psi.tilde.sim, phi.tilde.sim); invisible(gc())
rm(psi.tilde.sim); invisible(gc())


#############################################################################
# Top words --------------------------------------------------------------- #
#############################################################################

m = 15

#for Greek data only:
# lemmas.unique = read.csv("~/R files/Data/lemmas.unique.csv", header = TRUE, colClasses = "character")

# # top m words with highest posterior mean for each sense separately for each time period
# top.m.words = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k) {order(psi.tilde.post.mean[,k,t], decreasing = TRUE)[1:m]})}, simplify = "array")
# dimnames(top.m.words) = list(Rank = 1:m, Sense = 1:num.senses, Time = 1:num.periods)
# 
# top.m.LemmaIDs = aaply(top.m.words, 1:3, function(x) {words.used[x]})
# # noquote(aperm(top.m.LemmaIDs,c(2,1,3)))
# #for Greek data only:
# # top.m.original = aaply(top.m.LemmaIDs, 1:3, function(x) {lemmas.unique[match(x, lemmas.unique[,1]),2]})
# # noquote(aperm(top.m.original,c(2,1,3)))

# top m words with highest posterior mean for each sense across all time periods
psi.tilde.post.mean.marg.over.time = apply(psi.tilde.post.mean, 1:2, sum) / num.periods
top.m.words.marg.over.time = sapply(1:num.senses, function(k) {order(psi.tilde.post.mean.marg.over.time[,k], decreasing = TRUE)[1:m]})
dimnames(top.m.words.marg.over.time) = list(Rank = 1:m, Sense = 1:num.senses)

top.m.marg.over.time.LemmaIDs = aaply(top.m.words.marg.over.time, 1:2, function(x) {words.used[x]})
# noquote(t(top.m.marg.over.time.LemmaIDs))
#for Greek data only:
# top.m.marg.over.time.original = aaply(top.m.marg.over.time.LemmaIDs, 1:2, function(x) {lemmas.unique[match(x, lemmas.unique[,1]),2]})
# noquote(t(top.m.marg.over.time.original))

print("top words")
print(t(top.m.marg.over.time.LemmaIDs), quote = FALSE)
# print(t(top.m.marg.over.time.original), quote = FALSE) #for Greek data only


#############################################################################
# Normalised sense probabilities over true senses ------------------------- #
#############################################################################

#Specify number of true target-word senses
num.true.senses = 3

#Define sense order
#E.g. for kosmos, the 3 true senses are decoration, order, world
#We map model sense 2 to decoration; model sense 3 to order; and model senses 1 and 4 to world
#So we set desired.labelling = list(c(2), c(3), c(1,4))
desired.labelling = list(c(2), c(3), c(1,4))


#Normalise posterior sense probabilities over all senses
sense.probs.sim.normalised.all = sapply(1:num.snippets, function(d){
  sense.probs.sim[,,d]/rowSums(sense.probs.sim[,,d])}, simplify = "array")

#Posterior mean sense probabilities
sense.probs.sim.normalised = sense.probs.sim.normalised.all[,1:num.true.senses,]
for (k in 1:num.true.senses) {
  sense.probs.sim.normalised[,k,] = apply(sense.probs.sim.normalised.all[,desired.labelling[[k]],,drop=FALSE], 3, rowSums)
}
sense.probs.post.mean = apply(sense.probs.sim.normalised[(burn.in+1):num.samples,,], 2:3, mean)
rm(sense.probs.sim.normalised.all); invisible(gc())


#############################################################################
# Brier score ------------------------------------------------------------- #
#############################################################################

#Index of snippets to test
sense.ids = levels(factor(snippets.info$sense.id))
idx = which(snippets.info$sense.id %in% sense.ids[1:num.true.senses])
names(idx) = idx 
rm(sense.ids)

#Brier score
BS = multiclass.Brier(t(sense.probs.post.mean[,idx]), as.numeric(factor(snippets.info$sense.id[idx])))
print("brier score"); print(BS)


#############################################################################
# close for loop ---------------------------------------------------------- #
#############################################################################

}# for seed
