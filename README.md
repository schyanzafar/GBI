This repository contains the R scripts and data files used to produce the results reported in [Zafar and Nicholls (2024)](https://doi.org/10.48550/arXiv.2410.01475): "Exploring Learning Rate Selection in Generalised Bayesian Inference using Posterior Predictive Checks".

**Data:** This folder contains a separate subfolder for each of the target words used as a test case in the paper. Each subfolder contains the files _"snippets.RData"_ and _"embeddings.RData"_, which are required to run the EDiSC model on the relevant target word data. The English target word subfolders additionally contain the files _"snippets.words.RData"_, which show the actual English context words used in the snippets rather than the word IDs. The file _"lemmas.unique.csv"_ is a list of all the unique lemmas used in the Greek datasets, which is used to map the context word lemma IDs to the actual Greek words.

**Scripts:** This folder contains four files:
1. _"EDiSC_powerlik_model.stan"_ is the EDiSC model coded in Stan;
2. _"EDiSC_powerlik_model_setup.R"_ is the R script used to run the model;
3. _"Top_words_and_Brier_score.R"_ is the R script used to analyse the posterior for the highest probability context words under each model sense, and to compute the Brier score; and
4. _"PPC.R"_ is the R script used to carry out the posterior predictive checks, compute the p-values and produce the graphs.
