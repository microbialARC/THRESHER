library(TreeSim)
# Simulates phylogenetic trees using a birth-death process (TreeSim package).
#
# The birth-death tree represents the evolutionary history of a strain lineage,
# capturing diversification events regardless of whether they occur within a 
# single host or across multiple hosts through transmission. The process is 
# agnostic to host boundaries — it models the emergence and loss of sub-lineages
# over time.
#
# Parameters:
#   lambda (speciation rate): The rate at which new sub-lineages arise per 
#     existing lineage per unit time (year). Higher lambda produces more 
#     branching events and a bushier tree.
#
#   mu (extinction rate): The rate at which existing sub-lineages are lost
#     per lineage per unit time (year). This can reflect lineage extinction 
#     through bottlenecks, clearance, or competitive displacement.
#     When mu = 0, the model reduces to a pure-birth (Yule) process where 
#     all lineages persist — no lineage loss occurs.
#
# Note: lambda and mu have no default values and must be specified by the user
# in the configuration file.
# Function ----

phylogeny_simulation <- function(output_tree_path,
                                 simu_years = simu_years,
                                 simu_taxa = simu_taxa,
                                 lambda_rate,
                                 mu_rate,
                                 number_simu = 1,
                                 simu_seed){
  set.seed(simu_seed)
  
  simu_tree <- sim.bd.taxa.age(n=simu_taxa,
                               age = simu_years,
                               numbsim = number_simu,
                               lambda = lambda_rate,
                               mu = mu_rate,
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
lambda_rate <- as.numeric(snakemake@params[["lambda_rate"]])
mu_rate <- as.numeric(snakemake@params[["mu_rate"]])
simu_seed <- as.numeric(snakemake@params[["seed"]])
# Execute the function ----
# Each simulation only produces one tree, so number_simu is set to 1
phylogeny_simulation(output_tree_path = output_tree_path,
                     simu_years = simu_years,
                     simu_taxa = simu_taxa,
                     lambda_rate = lambda_rate,
                     mu_rate = mu_rate,
                     number_simu = 1,
                     simu_seed = simu_seed)
