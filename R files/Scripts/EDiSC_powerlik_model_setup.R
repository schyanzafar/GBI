options(digits = 4)
library(rstan)
rstan_options(auto_write = TRUE)

#Set working directory
setwd("")

#Read snippets data
snippets_dir = "~/R files/Data/apple"
snippets_name = "snippets.RData"
load(paste(snippets_dir, snippets_name, sep = "/"))

#Read embeddings matrix (rho)
embed_dir = "~/R files/Data/apple"
embed_name = "embeddings.RData"
load(paste(embed_dir, embed_name, sep = "/"))

#Specify model path
model_dir  = "~/R files/Scripts"
model_name = "EDiSC_powerlik_model.stan"
model_path = paste(model_dir, model_name, sep = "/")


#Set number of senses K
num.senses = 2

#Set coefficient alpha in the model X(t) = alpha * X(t-1) + noise
alpha.phi   = 0.9
alpha.theta = 0.9

#Set variances for phi, chi, theta and sigma
rho.squared.diffs = unlist(lapply(1:(nrow(rho)-1), function(x){as.numeric(colSums((rho[x,] - t(rho[-(1:x),,drop=FALSE]))^2))}))
median.rho.squared.diff = median(rho.squared.diffs)
kappa.phi   = 0.25
kappa.theta = 0.5 / median.rho.squared.diff
kappa.chi   = 2.5 / median.rho.squared.diff
kappa.sigma = 0.25
rm(rho.squared.diffs); gc()

#Set power to raise the likelihood to
lik.power = 1.0


#Prepare snippets to feed into Stan model
stan_snippets = snippets[,1:snippet.length] #discard extra columns
stan_snippets = t(apply(stan_snippets, 1, function(x){c(x[!is.na(x)], x[is.na(x)])})) #move NAs to end of each snippet
stan_snippets = replace(stan_snippets, is.na(stan_snippets), 0) #replace NAs with zero


#Prepare list of data to feed into Stan
stan_data = list(V = num.words,
                 K = num.senses,
                 G = num.genres,
                 T = num.periods,
                 D = num.snippets,
                 L = snippet.length, 
                 M = embedding.dim,
                 
                 snippets = stan_snippets, 
                 lengths  = snippet.lengths, 
                 times    = snippets$Time, 
                 genres   = snippets$genre,
                 
                 rho = rho,
                 
                 alpha_phi   = alpha.phi,
                 alpha_theta = alpha.theta, 
                 kappa_phi   = kappa.phi,
                 kappa_theta = kappa.theta, 
                 kappa_chi   = kappa.chi,
                 kappa_sigma = kappa.sigma,
                 
                 lambda = lik.power)

rm(stan_snippets); gc()


#Save data
save(stan_data, file = "EDiSC_fit_data.RData")
save(snippets, snippets.info, words.used, num.words, num.periods, num.genres,
     num.senses, num.snippets, snippet.length, snippet.lengths,  file = "snippets.RData")
save(rho, embedding.dim, file = "embeddings.RData")


#Fit model using Stan NUTS
seeds = c(100, 200)
target_accept_rate = 0.6
num_iterations = 2000
warmup = 1000
thin = 1

for (seed in seeds) {
  EDiSC_fit = stan(file = model_path, model_name = model_name, data = stan_data, verbose = FALSE,
                   pars = c("log_sense_probs", "phi_tilde", "psi_tilde", "theta", "chi", "sigma"),
                   control = list(metric = "diag_e", max_treedepth = 10, stepsize_jitter = 1,
                                  # adapt_init_buffer = 75, adapt_term_buffer = 50, adapt_window = 25, 
                                  adapt_delta = target_accept_rate), 
                   cores = 1, chains = 1, iter = num_iterations, warmup = warmup, thin = thin, seed = seed)
  
  #Save results
  save(EDiSC_fit, file = paste0("EDiSC_fit_seed_",seed,".RData"))
}

