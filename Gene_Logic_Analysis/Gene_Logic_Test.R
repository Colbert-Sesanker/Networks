# Main Logic Analysis Function
# Colbert Sesanker 2013

source("Gene_Logic.R")

default_triplet = c("Endo16_1","IrxA_1","gsc_1")
logic_analysis <- function (triplet=default_triplet, sig_measure="R2", sample_density=3, time=1, logic="OR") { 
    results <- test_logic(triplet[1],triplet[2],triplet[3], sig_measure, sample_density)
    visualize_logic(results,time,logic)
    dev.new()
    logic_hist(results,time,logic)
    dev.new()
    scatter_logic(results,time,logic)
}

# Tests random sources to a given target
genes   <- featureNames(GeneData)
rand_analysis <- function(target="gsc_1", time=1) {
    triplet <- c(sample(genes,1), sample(genes,1), target, time=time) 
    logic_analysis(triplet)
}
