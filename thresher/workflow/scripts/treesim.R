library(TreeSim)
# Use the parameters in the manuscript for Epidemic simulations
# https://www.sciencedirect.com/science/article/pii/S1755436524000562
# In the epidemiological context
# lambda (birth rate) is the rate at which an infected individual transmits the disease to an uninfected individual 
# mu (death rate) is the rate at which an infected individual becomes non-infectious

# Use 
# birth_rate = 1.5
# death_rate = 1
# Function ----
phylogeny_simulation <- function(output_tree_path,
                                 simu_years = simu_years,
                                 simu_taxa = simu_taxa,
                                 birth_rate = 1.5,
                                 death_rate = 1,
                                 number_simu = 1,
                                 simu_seed){
  set.seed(simu_seed)
  
  simu_tree <- sim.bd.taxa.age(n=simu_taxa,
                               age = simu_years,
                               numbsim = number_simu,
                               lambda = birth_rate,
                               mu = death_rate,
                               frac = 1,
                               mrca = TRUE)
  
  # export the tree 
  write.tree(simu_tree,
             file = output_tree_path)
  
}
# Import from Snakemake ----
output_tree_path <- snakemake@output[["sim_tree"]]
simu_years <- as.numeric(snakemake@params[["years"]])
simu_taxa <- as.numeric(snakemake@params[["taxa"]])
simu_seed <- as.numeric(snakemake@params[["seed"]])
# Execute the function ----
phylogeny_simulation(output_tree_path = output_tree_path,
                     simu_years = simu_years,
                     simu_taxa = simu_taxa,
                     birth_rate = 1.5,
                     death_rate = 1,
                     number_simu = 1,
                     simu_seed = simu_seed)
