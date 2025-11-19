import math
import re
import numpy
import pandas as pd
import sys
import random
import os
import glob
from Bio import SeqIO
import subprocess

# The SequenceMutator class encapsulates logic for simulating evolution.
# It supports both uniform random mutation and weighted mutation based on entropy data.
# - If use_weight_mutate is True and a valid CSV is provided, mutations are sampled according to
#   normalized entropy values for genome bins, allowing for region-specific mutation rates.
# - Otherwise, mutations are applied uniformly at random across the sequence.
# - The class uses a probability matrix derived from the chosen substitution model (e.g., JC69, K2P, K3P, GTR)
#   to determine possible base substitutions.
# - Each mutation records detailed information (node, position, original and mutated base) for downstream analysis.
class SequenceMutator:
	def __init__(self, probability, use_weight_mutate, weight_mutate_data, concatenated_contig):
		self.probability = probability
		self.use_weight_mutate = use_weight_mutate
		self.weight_mutate_data = weight_mutate_data
		self.concatenated_contig = concatenated_contig
		# Read weighted mutation data in csv format if use_weight_mutate is True
		if self.use_weight_mutate and self.weight_mutate_data is not None:
			try:
				self.weight_mutation_df = pd.read_csv(self.weight_mutate_data)
				# Ensure the DataFrame is in the correct format
				if not all(col in self.weight_mutation_df.columns for col in ['position', 'entropy']):
					raise ValueError("Weighted mutation data is missing required columns. Required: position, entropy")
				# Select only position and entropy columns
				self.weight_mutation_df = self.weight_mutation_df[['position', 'entropy']]
				# There are 2 types of positions here. The first is the alignment position while the second is the original entropy profile position
				# The alignment position is the position in the MSA, while the profile position is the position in the original start genome for simulation
				# The alignment position is mutatable due to gain and loss of MGEs, while the profile position is fixed and will be used for visualization
				self.weight_mutation_df['profile_position'] = self.weight_mutation_df['position'].astype(int)
				self.weight_mutation_df['alignment_position'] = self.weight_mutation_df['position'].astype(int)
				# Remove the original 'position' column after assigning alignment and profile positions
				self.weight_mutation_df.drop(columns=['position'], inplace=True)
				self.entropy_value = self.weight_mutation_df['entropy'].copy()
				print("Use weighted mutation.")
				# Check if entropy value of 0 is found, if so the corrected minimal basal value will be applied
				if self.entropy_value.min() > 0:
					print(f"Minimal entropy value found: {self.entropy_value.min():.4e}")
					print(f"No zero entropy values found. No correction will be applied.")
				else:
					# For the entropy_value <= 0, use the minimal entropy value as basal probability
					# This ensures that every position has a non-zero chance of mutation based on the hypothesis that every site can potentially mutate
					# This might not be absolutely biologically accurate, but it prevents positions from being completely immutable in the simulation
					self.minimal_entropy = self.entropy_value[self.entropy_value > 0].min()
					self.entropy_value[self.entropy_value <= 0] = self.minimal_entropy
					print(f"Zero entropy values found. Correction will be applied.")
					print(f"Corrected minimal basal value for zero entropy: {self.entropy_value.min():.4e}")
				# Normalize entropy values to probabilities using linear scaling and add back to weight_mutation_df
				print("Normalising entropy values.")
				self.normalized_entropy_value = self.entropy_value / self.entropy_value.sum()
				self.weight_mutation_df['normalized_entropy'] = self.normalized_entropy_value
				
			# If fail reading the csv, rolling back to uniform random mutation
			except Exception as e:
				print(f"Error reading weighted mutation data: {e}")
				print("Rolling back to uniform random selection")
				self.weight_mutation_df = None
				use_weight_mutate = False
		else:
			# If weighted mutation is not used, create a data frame to map the transient alignment position to profile position mapping
			# use concatenated_contig as the reference sequence
			self.position_df = pd.DataFrame({
				'alignment_position': range(len(self.concatenated_contig)),
				'profile_position': range(len(self.concatenated_contig))
			})

	# This is the function for weighted mutation using the entropy values
	def weight_mutate(self, seq, node = None):
		print(f"Simulating weighted mutation for {node}.")
		# Convert to list for in-place mutation
		out = list(seq)
		# Weighted random selection using normalized_entropy_value
		if self.use_weight_mutate and self.weight_mutate_data is not None:
			# Select a random alignment position using normalized entropy values as weights
			pos = numpy.random.choice(self.weight_mutation_df['alignment_position'], p=self.weight_mutation_df['normalized_entropy'])
			original_base = out[pos]
			# If the original base is '-', meaning that this location is lost in the genome
			# Keep selecting new positions until the original base is not '-'
			while original_base == "-":
				pos = numpy.random.choice(self.weight_mutation_df['alignment_position'], p=self.weight_mutation_df['normalized_entropy'])
				original_base = out[pos]
			# probability is determined by the substitution model
			new_base = numpy.random.choice(self.probability[original_base])
			# instead of slicing and concatenating strings
			# directly modify in place in the list
			out[pos] = new_base
			# Store mutation information for retrieval
			# Profile position is used to generate the visualization
			# The alignment position only indicates the transient position in the MSA at the time of mutation
			# Any subsequent gain or loss of MGEs may shift the alignment position
			mutation_info = {
				'node': str(node),
				'alignment_position': int(pos),
				'profile_position': int(self.weight_mutation_df[self.weight_mutation_df['alignment_position'] == pos]['profile_position'].values[0]),
				'original_base': str(original_base),
				'mutated_base': str(new_base)
			}
			# return both mutated sequence and info at once
			return ''.join(out), mutation_info
		else:
			# If weight mutation is not available, fall back to uniform mutation
			print("WARNING: Weighted mutation is not available. Falling back to uniform mutation.")
			return self.uniform_mutate(seq, node)

	# This is the function using uniform random mutation 
	def uniform_mutate(self, seq, node = None):
		print(f"Simulating uniform mutation for {node}.")
		seq_len = len(seq)
		# Convert to list for in-place mutation
		out = list(seq)
		# this section is the range where the mutation could happen in the genome
		selection = range(seq_len)
		# Within the count of mutation, randomly set mutations from the genome
		pos = numpy.random.choice(selection)
		original_base = out[pos]
		# Same reason for while loop mentioned above
		while original_base == "-":
			pos = numpy.random.choice(selection)
			original_base = out[pos]
		# probability is determined by the substitution model
		new_base = numpy.random.choice(self.probability[original_base])
		# instead of slicing and concatenating strings
		# which could be slow
		# I just modify in place
		out[pos] = new_base
		# Store mutation information for retrieval
		# Get the profile position corresponding to the selected alignment position
		# Get the profile position from position_df if weight mutation is not used
		mutation_info = {
			'node': str(node),
			'alignment_position': int(pos),
			'profile_position': int(self.position_df[self.position_df['alignment_position'] == pos]['profile_position'].values[0]),
            'original_base': str(original_base),
            'mutated_base': str(new_base)
        }

		# return both mutated sequence and info at once
		return ''.join(out), mutation_info

# SequenceMGEser simulates the gain and loss of MGEs in the genome
class SequenceMGEser:
	def __init__(self, chr_bins, gain_rate, loss_rate, mge_fasta, mge_entropy):
		self.gain_rate = gain_rate
		self.loss_rate = loss_rate
		self.chr_bins = chr_bins
		self.mge_fasta = mge_fasta
		self.mge_entropy = mge_entropy
	# Read the chr_bins to get the locations of all bins and whether or not they are MGEs
	# If they are not MGE, the genomes are NOT allowed to lose the bins
		# Read the bin data
		if self.gain_rate and self.loss_rate and self.chr_bins is not None:
			try:
				self.chr_bins_df = pd.read_csv(self.chr_bins)
				self.chr_bins_df["bin_index"] = self.chr_bins_df["bin_index"].astype(str)
				# Error if wrong format
				if not all(col in self.chr_bins_df.columns for col in ['bin_index', 'gene', 'start', 'end', 'mge_id']):
					raise ValueError("Bin data is missing required columns. Required: bin_index, gene, start, end, and mge_id")
				# Select only needed columns
				self.chr_bins_df = self.chr_bins_df[['bin_index', 'gene', 'start', 'end', 'mge_id']]
				# Also there are 2 types of start and end postions here (Alignment and Profile)
				# The alignment start and end positions are the positions in the MSA, while the profile start and end positions are the positions in the original genome
				# The alignment start and end positions are mutatable due to gain and loss of MGEs, while the profile start and end positions are fixed and will be used for visualization
				self.chr_bins_df['alignment_start'] = self.chr_bins_df['start'].astype(int)
				self.chr_bins_df['alignment_end'] = self.chr_bins_df['end'].astype(int)
				self.chr_bins_df['profile_start'] = self.chr_bins_df['start'].astype(int)
				self.chr_bins_df['profile_end'] = self.chr_bins_df['end'].astype(int)
				# Remove the original 'start' and 'end' columns after assigning alignment and profile positions
				self.chr_bins_df.drop(columns=['start', 'end'], inplace=True)
				# The bins within the chromosome which allow gene loss
				self.loss_bin = self.chr_bins_df['bin_index'][self.chr_bins_df['mge_id'] != "non-mge"].astype(str).tolist()
				# The bins within the chromosome which allow insertion (gene Gain)
				self.insertion_bin = self.chr_bins_df['bin_index'][self.chr_bins_df['gene'] == 'non-cds'].astype(str).tolist()
			except Exception as e:
				print(f"MGE gain and loss simulation is not available due to error reading bin data: {e}")
		else:
			print("Provided bin data is not complete.")
		# Read the MGE database to get the sequences of all MGEs for gene gain simulation
		if self.mge_fasta and self.mge_entropy is not None:
			try:
				self.mge_fasta_dict = SeqIO.to_dict(SeqIO.parse(self.mge_fasta, "fasta"))
				self.mge_fasta_dict_idx = list(self.mge_fasta_dict.keys())
				# search for all entropy files in self.mge_entropy and get a dict with keys as mge ids and values as entropy values
				self.mge_entropy_dict = {
					re.sub(r"_entropy\.csv$", "", os.path.basename(entropy_entry)): entropy_entry
					for entropy_entry in glob.glob(os.path.join(self.mge_entropy, "*_entropy.csv"))
				}
			except Exception as e:
				print(f"Error loading MGE database: {e}")

	def mge_gain(self, seq, node = None):
		# For gene gain simulation
		print(f"Simulating MGE gain for {node}.")
		# Safety check
		if not hasattr(self, 'insertion_bin') or not self.insertion_bin:
			print(f"Error: No insertion bins available for node {node}")
			return seq, None
		
		if len(self.insertion_bin) == 0:
			print(f"Error: No insertion bins available for node {node}")
			return seq, None
		try:
			# First the bin containing the insertion is randomly chosen
			bin_gain_index = str(numpy.random.choice(self.insertion_bin))
		except Exception as e:
			print(f"Error choosing insertion bin for node {node}: {e}")
			print(f"insertion_bin state: {self.insertion_bin}")
			return seq, None
		
		bin_gain_row = self.chr_bins_df.index[self.chr_bins_df['bin_index'] == bin_gain_index][0]
		bin_gain_mge_id = str(self.chr_bins_df['mge_id'][self.chr_bins_df['bin_index'] == bin_gain_index].values[0])
		# The position of insertion is also randomly picked within the bin
		bin_gain_aln_start = self.chr_bins_df[self.chr_bins_df["bin_index"] == bin_gain_index]["alignment_start"].values[0]
		bin_gain_aln_end = self.chr_bins_df[self.chr_bins_df["bin_index"] == bin_gain_index]["alignment_end"].values[0]
		bin_gain_aln_length = bin_gain_aln_end - bin_gain_aln_start + 1

		# if the length of the insertion bin is less than 2, keep picking new bins until a valid one is found
		# because a bin with length less than 2 cannot accommodate an insertion
		while bin_gain_aln_length <= 2:
			try:
				bin_gain_index = str(numpy.random.choice(self.insertion_bin))
			except Exception as e:
				print(f"Error choosing insertion bin for node {node}: {e}")
				print(f"insertion_bin state: {self.insertion_bin}")
				return seq, None
			bin_gain_row = self.chr_bins_df.index[self.chr_bins_df['bin_index'] == bin_gain_index][0]
			bin_gain_mge_id = str(self.chr_bins_df['mge_id'][self.chr_bins_df['bin_index'] == bin_gain_index].values[0])
			# The position of insertion is also randomly picked within the bin
			bin_gain_aln_start = self.chr_bins_df[self.chr_bins_df["bin_index"] == bin_gain_index]["alignment_start"].values[0]
			bin_gain_aln_end = self.chr_bins_df[self.chr_bins_df["bin_index"] == bin_gain_index]["alignment_end"].values[0]
			bin_gain_aln_length = bin_gain_aln_end - bin_gain_aln_start + 1
		try:
			# 2 insertion posstions: alignment and profile(for visualization)
			bin_gain_aln_insertion = numpy.random.choice(range(bin_gain_aln_start, bin_gain_aln_end))
			# calculate the distance of insertion position to the start of the bin
			dist_start = bin_gain_aln_insertion - bin_gain_aln_start
			# use the distance to calculate the profile insertion position
			bin_gain_profile_insertion = self.chr_bins_df[self.chr_bins_df["bin_index"] == bin_gain_index]["profile_start"].values[0] + dist_start
		except Exception as e:
			print(f"Error choosing insertion position within bin {bin_gain_index} for node {node}: {e}")
			print(f"bin_gain_aln_start: {bin_gain_aln_start}, bin_gain_aln_end: {bin_gain_aln_end}")
			return seq, None
		
		try:
			# The MGE to be inserted is randomly picked from the MGE database
			insert_mge_index = numpy.random.choice(list(self.mge_fasta_dict_idx))
		except Exception as e:
			print(f"Error choosing MGE to insert for node {node}: {e}")
			return seq, None
		# Get the sequence of the MGE to be inserted and save it as string
		insert_mge_seq = ''.join(self.mge_fasta_dict[insert_mge_index].seq)
		# The length of the to be insert MGE
		length_insert_mge_seq = len(insert_mge_seq)
		# Get the entropy values for this MGE
		insert_mge_entropy = pd.read_csv(self.mge_entropy_dict[insert_mge_index])
		# Update the new sequence
		# If the sequence is a list, join it to string first
		if isinstance(seq, list):
			seq = ''.join(seq)
		new_seq = seq[:bin_gain_aln_insertion] + insert_mge_seq + seq[bin_gain_aln_insertion:]

		# Also update all existing sequences, because the exported fasta file containing the MSA of all simulated sequences
		# And the entropy values are applied to every sequence
		for key_entry in sequence.keys():
			if key_entry != str(node):
				# Mark the inserted region with '-' to keep the sequence length the same
				seq_updated =  sequence[key_entry][:bin_gain_aln_insertion] + "-" * length_insert_mge_seq + sequence[key_entry][bin_gain_aln_insertion:]
				sequence[key_entry] = seq_updated

		# Update chromosome bin positions after MGE insertion.
		# Every time an MGE is inserted, the alignment positions of all bins after the insertion site need to be shifted by the length of the inserted MGE.
		# Rows above bin_gain_row remain unchanged.
		# The original insertion bin row is excluded, since it is replaced by new bins.
		# Lower bins (after the insertion) have their start and end positions shifted to account for the inserted MGE length.

		upper_chr_bins_df = self.chr_bins_df.iloc[:bin_gain_row]
		lower_chr_bins_df = self.chr_bins_df.iloc[bin_gain_row+1:].copy()
		lower_chr_bins_df["alignment_start"] += length_insert_mge_seq
		lower_chr_bins_df["alignment_end"] += length_insert_mge_seq

		# there are 3 rows for insert_chr_bins_df
		# 2 rows are the inserted bin being split into 2 parts
		# 1 row is the new inserted MGE
		insert_chr_bins_df = pd.DataFrame({
			"bin_index": [f"{bin_gain_index}_insert_1", f"{bin_gain_index}_insert_2", f"{bin_gain_index}_insert_3"],
			"gene": [f"{self.chr_bins_df.iloc[bin_gain_row]['gene']}_1", f"{insert_mge_index}", f"{self.chr_bins_df.iloc[bin_gain_row]['gene']}_2"],
			"alignment_start": [bin_gain_aln_start, 
					   bin_gain_aln_insertion + 1,
					   bin_gain_aln_insertion + length_insert_mge_seq + 1],
			"alignment_end": [bin_gain_aln_insertion,
					 bin_gain_aln_insertion + length_insert_mge_seq,
					 bin_gain_aln_end + length_insert_mge_seq],
			# For profile start and end.
			# The bin containing the insertion is split into 2 parts
			# The profile positions of the original bin remain unchanged
			# All sites in the inserted MGE takes the same profile position at the insertion site

			
			"profile_start": [self.chr_bins_df.iloc[bin_gain_row]['profile_start'],
					 bin_gain_profile_insertion,
					 bin_gain_profile_insertion + 1],
			"profile_end": [bin_gain_profile_insertion - 1,
				   bin_gain_profile_insertion,
				   self.chr_bins_df.iloc[bin_gain_row]['profile_end']],
			"mge_id": [bin_gain_mge_id,insert_mge_index,bin_gain_mge_id]
		})

		# Update the loss_bin to add _insert_2
		self.loss_bin.append(f"{bin_gain_index}_insert_2")
		# If the bin containing the insertion is an MGE, then _insert_1 and _insert_3 can also be lost
		if bin_gain_mge_id != "non-mge":
			self.loss_bin.append(f"{bin_gain_index}_insert_1")
			self.loss_bin.append(f"{bin_gain_index}_insert_3")
		# Update the insertion_bin to add _insert_1, _insert_2, _insert_3 and remove the original bin index
		self.insertion_bin.append(f"{bin_gain_index}_insert_1")
		self.insertion_bin.append(f"{bin_gain_index}_insert_2")
		self.insertion_bin.append(f"{bin_gain_index}_insert_3")
		
		# If the bin is in insertion_bin, remove the original bin index
		if bin_gain_index in self.insertion_bin:
			self.insertion_bin.remove(bin_gain_index)
		if bin_gain_index in self.loss_bin:
			self.loss_bin.remove(bin_gain_index)
		
		# Update the chr_bins_df
		self.chr_bins_df = pd.concat([upper_chr_bins_df, insert_chr_bins_df, lower_chr_bins_df], ignore_index=True)

		# Update the weight_mutation_df in mutator
		print("Updating weighted mutation data after MGE insertion.")
		# upper part remains the same
		upper_part = mutator.weight_mutation_df[['alignment_position','profile_position', 'entropy']][mutator.weight_mutation_df['alignment_position'] < bin_gain_aln_insertion].copy()
		# alignment positions in the lower part are shifted by length_insert_mge_seq
		# profile positions in the lower part remain unchanged
		lower_part = mutator.weight_mutation_df[['alignment_position','profile_position', 'entropy']][mutator.weight_mutation_df['alignment_position'] >= bin_gain_aln_insertion].copy()
		lower_part['alignment_position'] += length_insert_mge_seq
		# new part is the entropy values of the inserted MGE
		new_part = pd.DataFrame({
			'alignment_position': range(bin_gain_aln_insertion, bin_gain_aln_insertion + length_insert_mge_seq),
			# profile positions in the new part are set to the insertion profile position
			'profile_position': bin_gain_profile_insertion,
			'entropy': insert_mge_entropy['entropy'].values
		})
		# Concatenate all parts
		updated_weight_mutation_df = pd.concat([upper_part, new_part, lower_part], ignore_index=True)
		# Re-order by alignment_position
		updated_weight_mutation_df = updated_weight_mutation_df.sort_values(by='alignment_position').reset_index(drop=True)
		# Re-normalize the entropy values
		updated_entropy_value = updated_weight_mutation_df['entropy'].copy()
		# Check if entropy value <= 0 is found, if so the corrected minimal basal value will be applied
		if updated_entropy_value.min() > 0:
			print(f"Minimal entropy value found after MGE insertion: {updated_entropy_value.min():.4e}")
			print(f"No zero entropy values found. No correction will be applied.")
		else:
			# For the entropy_value <= 0, add basal probability using the minimal entropy value
			minimal_entropy = updated_entropy_value[updated_entropy_value > 0].min()
			updated_entropy_value[updated_entropy_value <= 0] = minimal_entropy
			print(f"Zero entropy values found after MGE insertion. Correction will be applied.")
			print(f"Corrected minimal basal value for zero entropy: {updated_entropy_value.min():.4e}")
		print("Re-normalising entropy values after MGE insertion.")
		updated_normalized_entropy_value = updated_entropy_value / updated_entropy_value.sum()
		updated_weight_mutation_df['normalized_entropy'] = updated_normalized_entropy_value

		# Update the weight_mutation_df in mutator
		mutator.weight_mutation_df = updated_weight_mutation_df
		# Store gene gain information for retrieval
		gene_gain_info = {
			'node': str(node),
			'insert_bin': str(bin_gain_index),
			'insertion_profile_position': int(bin_gain_profile_insertion),
			'insert_mge_id': str(insert_mge_index),
			'insert_mge_length': int(length_insert_mge_seq)
		}
		# return both mutated sequence and info at once
		return new_seq, gene_gain_info
	
	# Function for MGE loss simulation
	def mge_loss(self, seq, node = None):
		# For gene loss simulation
		print(f"Simulating MGE loss for {node}.")
		# Perform MGE loss simulation steps here
		# First the bin to be lost is randomly chosen from self.loss_bin
		bin_loss_index = str(numpy.random.choice(self.loss_bin))
		bin_loss_mge_id = str(self.chr_bins_df['mge_id'][self.chr_bins_df['bin_index'] == bin_loss_index].values[0])
		# The position of loss is the whole bin
		bin_loss_start = self.chr_bins_df[self.chr_bins_df["bin_index"] == bin_loss_index]["alignment_start"].values[0]
		bin_loss_end = self.chr_bins_df[self.chr_bins_df["bin_index"] == bin_loss_index]["alignment_end"].values[0]
		# Extract the bin sequence to check if there is already a loss (all '-')
		bin_loss_seq = seq[bin_loss_start:bin_loss_end]
		# If the bin sequence is all '-', meaning that this location is already lost in the
		# genome, keep selecting new bins while the bin sequence is all '-'
		while all(base == '-' for base in bin_loss_seq):
			bin_loss_index = str(numpy.random.choice(self.loss_bin))
			bin_loss_mge_id = str(self.chr_bins_df['mge_id'][self.chr_bins_df['bin_index'] == bin_loss_index].values[0])
			bin_loss_start = self.chr_bins_df[self.chr_bins_df["bin_index"] == bin_loss_index]["alignment_start"].values[0]
			bin_loss_end = self.chr_bins_df[self.chr_bins_df["bin_index"] == bin_loss_index]["alignment_end"].values[0]
			bin_loss_seq = seq[bin_loss_start:bin_loss_end]
		# The length of the to be lost MGE
		length_loss_mge_seq = bin_loss_end - bin_loss_start
		# Update the new sequence
		# Mark the lost region with '-' to keep the sequence length the same
		new_seq = seq[:bin_loss_start] + "-" * (bin_loss_end - bin_loss_start) + seq[bin_loss_end:]
		# Store gene loss information for retrieval
		gene_loss_info = {
			'node': str(node),
			'loss_bin': str(bin_loss_index),
			'loss_mge_id': str(bin_loss_mge_id),
			'loss_mge_length': int(length_loss_mge_seq)
		}
		return new_seq, gene_loss_info

################# Analysis parameters imported from snakemake#################

# Input 
prefix = str(snakemake.params.prefix)                        		# Prefix of the output files
output_path = str(snakemake.params.output_dir)      				# Path to the output directory
intermediate_path = str(snakemake.params.intermediate_dir)      	# Path to the intermediate files directory
path_to_seq = str(snakemake.input.ancestor)      		    		# Path to the starting ancestor genome

# Substitution Simulation Parameters
sub_model = str(snakemake.params.substitution_model)      			# Substitution Model
sub_rate = float(snakemake.params.mutation_rate)  			        # Substitution Rate (Substitutions/site/year)
sub_param = str(snakemake.params.model_parameters)        	        # Parameters of substitution model
use_weight_mutate = bool(snakemake.params.use_weighted_mutation)  	# Whether to use weighted mutation
weight_mutate_data = snakemake.params.weighted_mutation_data 		# Weighted mutation data

# Recombination Simulation Parameters
use_recomb = bool(snakemake.params.use_recombination) 				# Whether to use recombination in simulation
recom_rate = float(snakemake.params.recombination_rate) 			# Recombination rate (Recombination events/year)
mean_recomb_size = float(snakemake.params.mean_recombination_size) 	# Mean recombination size (bp)
min_recomb_size = float(snakemake.params.min_recombination_size)	# Minimal recombination size (bp)
nu = float(snakemake.params.nu)                						# nu in recombination simulation (snps/total length of recombinant fragment)

# MGEs Gain and Loss Simulation Parameters
use_mge = bool(snakemake.params.use_gain_loss) 						# Whether to use MGE in simulation
gain_rate = float(snakemake.params.gain_rate) 						# Gain rate (Gain events/year)
loss_rate = float(snakemake.params.loss_rate) 						# Loss rate (Loss events/year)
chr_bins = str(snakemake.params.chromosomal_bins)     				# Chromosomal bins file
mge_fasta = str(snakemake.params.mge_fasta)       					# Fasta file containing all MGE sequences
mge_entropy = str(snakemake.params.mge_entropy)       				# The directory containing entropy of each MGE
# Reproducibility
rseed = int(snakemake.params.seed)             						# Seed

################# Perform Analysis #################

# Set the seed
random.seed(rseed)

alpha=["A","T","C","G"]
transition={}
transition["A"]="G"
transition["G"]="A"
transition["C"]="T"
transition["T"]="C"
other1,other2={},{}
other1["A"],other2["A"]="C","T"
other1["T"],other2["T"]="G","A"
other1["C"],other2["C"]="G","A"
other1["G"],other2["G"]="C","T"


# Creates the substitution model with Kappa or sub_param
probability={}
if sub_model == "JC69":
	for N in alpha:
		toto = list(alpha)
		toto.remove(N)
		probability[N] = toto
elif sub_model == "K2P":
	# If the model is K2P
	# sub_param is Kappa
	kappa = float(sub_param)
	for N in alpha:
		probability[N]=[]
		i=0
		while i < round(10 * kappa,0):
			probability[N].append(transition[N])
			i+=1
		i=0
		while i < 5:
			probability[N].append(other1[N])
			probability[N].append(other2[N])
			i+=1
elif sub_model == "K3P":
	beta={}
	beta["A"]="C"
	beta["C"]="A"
	beta["G"]="T"
	beta["T"]="G"
	gamma={}
	gamma["A"]="T"
	gamma["T"]="A"
	gamma["C"]="G"
	gamma["G"]="C"
	toto = sub_param.split(",")
	ts,tv1,tv2= float(toto[0].strip(" ")) ,  float(toto[1].strip(" "))  , float(toto[2].strip(" ")) 
	for N in alpha:
		probability[N]=[]
		i=0
		while i < round(100 * ts,0):
			probability[N].append(transition[N])
			i+=1
		i=0
		while i < round(100 * tv1,0):
			probability[N].append(beta[N])
			i+=1
		i=0
		while i < round(100 * tv2,0):
			probability[N].append(gamma[N])
			i+=1
elif sub_model == "GTR":
	toto = sub_param.split(",")
	a,b,c,d,e,f = float(toto[0].strip(" ")) ,  float(toto[1].strip(" "))  , float(toto[2].strip(" ")) ,float(toto[3].strip(" ")) ,  float(toto[4].strip(" "))  , float(toto[5].strip(" ")) 
	for N in alpha:
		probability[N]=[]
	i=0
	while i < round(100 * a,0):
		probability["A"].append("G")
		probability["G"].append("A")
		i+=1
	i=0
	while i < round(100 * b,0):
		probability["A"].append("C")
		probability["C"].append("A")
		i+=1
	i=0
	while i < round(100 * c,0):
		probability["A"].append("T")
		probability["T"].append("A")
		i+=1
	while i < round(100 * d,0):
		probability["G"].append("C")
		probability["C"].append("G")
		i+=1
	while i < round(100 * e,0):
		probability["G"].append("T")
		probability["T"].append("G")
		i+=1
	while i < round(100 * f,0):
		probability["C"].append("T")
		probability["T"].append("C")
		i+=1

# In my version
# The function to generate random sequence is removed
# The starting genome must be provided
# ncontigs stores the number of contigs in the starting genome
ncontig=0
# seq stores the characters of starting genome
concatenated_contig = ""
with open(path_to_seq, "r") as seq_file:
	# Count the number of sequences
	file_content = seq_file.read()
	header_count = file_content.count(">")

	if use_weight_mutate and header_count > 1:
		# Weighted mutation simulation cannot be used with a genome containing multiple contigs.
		use_weight_mutate = False
		print("WARNING: Weighted mutation can not be used with a multi-contig genome and has been disabled.")

	if header_count > 1:
		print(f"Genome with {header_count} contigs provided, which will be concatenated")

		lines = file_content.strip().split('\n')

		for line in lines:
			# Remove all whitespace (spaces, tabs, newlines)
			line = line.strip()
			if not line:
				# Skip empty lines after removing all whitespace
				continue 

			if line.startswith(">"): # Header line detected
				ncontig += 1
				print(f"Concatenating #{ncontig} contig in the genome: {line}")
			else:
				# clean and validate line
				clean_line = line.upper() 
				# Remove any remaining whitespace and validate nucleotides
				clean_sequence = re.sub(r'[^ATCG]', '', line.upper())
				concatenated_contig += clean_sequence

    	# Export cleaned sequence as FASTA file
		with open(intermediate_path + "cleaned_sequence.fasta", "w") as output_file:
			output_file.write(">Cleaned_sequence\n")
			# Write sequence in lines of 80 characters (standard FASTA format)
			for i in range(0, len(concatenated_contig), 80):
				output_file.write(concatenated_contig[i:i+80] + "\n")
			print("Cleaned sequence exported to output directory as cleaned_sequence.fasta")
	else:
		print("The provided genome is a single-contig genome.\nNo cleaned_sequence.fasta will be generated.")
		lines = file_content.strip().split('\n')
		for line in lines:
			line = line.strip()
			if not line or line.startswith('>'):
				continue
			# Use regex for nucleotide parsing
			clean_sequence = re.sub(r'[^ATCG]', '', line.upper())
			concatenated_contig += clean_sequence
			# Only perform concatenation and export cleaned_sequence.fasta if not using weighted mutation
			if not use_weight_mutate and ncontig > 1:
				with open(intermediate_path + "cleaned_sequence.fasta", "w") as output_file:
					output_file.write(">Cleaned_sequence\n")
					for i in range(0, len(concatenated_contig), 80):
						output_file.write(concatenated_contig[i:i+80] + "\n")
				print("Cleaned sequence exported to output directory as cleaned_sequence.fasta")

# Activate the mutator function with probability matrix
mutator = SequenceMutator(
	probability=probability,
	use_weight_mutate=use_weight_mutate,
	weight_mutate_data=weight_mutate_data,
	concatenated_contig=concatenated_contig
)

mgeser = SequenceMGEser(
	chr_bins=chr_bins,
	gain_rate=gain_rate,
	loss_rate=loss_rate,
	mge_fasta=mge_fasta,
	mge_entropy=mge_entropy
)

branch_type = {}
branch_length = {}
with open(intermediate_path +"renamed.txt","r") as f:
	for l in f:
		a=l.strip("\n").split("\t")
		node = a[0]
		length = float(a[2])
		branch_length[node] = length
		branch_type[node] = a[3]


parent={}
dichotomies={}
with open(intermediate_path + "dichotomies.txt","r") as f:
	for l in f:
		a=l.strip("\n").split("\t")
		dichotomies[a[0]] = [a[1],a[2]]
		parent[a[1]] = a[0]
		parent[a[2]] = a[0]

roots=[]
with open(intermediate_path + "roots.txt","r") as f:
	for l in f:
		roots.append(l.strip("\n"))

##### Analysis of the time-scale phylogeny #####
tmp = list(roots)
cumul={}

for node in roots:
	cumul[node] = branch_length[node]

# Remove the 1==1 and some temporary codes
# this code traverses a phylogenetic tree to simulate sequence evolution from root to tips
# for each internal node, generates mutated sequences for child nodes based on branch lengths
# and accumulates evolutionary distances from root to each node

while len(tmp)>0:
	node = tmp[0]
	# process internal nodes by creating mutated sequences for their children
	if branch_type[node] != "tip":
		node1,node2 = dichotomies[node][0],dichotomies[node][1]
		tmp.append(node1)
		tmp.append(node2)
		# track cumulative evolutionary distance from root to each node
		if node in cumul:
			cumul[node1] = cumul[node] + branch_length[node1]
			cumul[node2] = cumul[node] + branch_length[node2]
		else:
			cumul[node1] =  branch_length[node1]
			cumul[node2] =  branch_length[node2]
		tmp.remove(node)
	else:
		# tip nodes are terminal. no further processing needed
		# new = mutate(sequence[node] , mutation_rate, branch[node])
		# remove new in here
		# this new doesn't seem to go anywhere 
		tmp.remove(node)

# Group nodes by evolutionary distance
rev={}
for node in cumul:
	if cumul[node] in rev:
		rev[cumul[node]].append(node)
	else:
		rev[cumul[node]]=[node]

values = list(rev.keys())
values.sort()

# Build the dictionary interval to define which branches overlapped in time
interval={}
for node in cumul:
	interval[node] = [cumul[node] - branch_length[node] , cumul[node] ]

# Creat new internal nodes
# The original code set MIN,MAX=10,0
# In my version, the unit of the branch length of the tree is year 
# That say, if any branch length is ≥ 10
# the MIN won't be updated correctly.
# My correction:
# the min_period is set to max in the system and will be corrected
# the max_period is set to 0 and will be corrected
min_period = sys.maxsize
max_period = 0
for node in cumul:
	if cumul[node]>max_period:
		max_period = cumul[node]
	if branch_length[node] != 0:
		if branch_length[node] < min_period:
			min_period= branch_length[node]


# Time-slice scanning to find overlapping branches
# This is the input for determining the time intervals downstream
branch_overlap_windows = {}
current_interval = 0
time_cursor = 0

while time_cursor < max_period:
	tmp=[]
	for node in interval:
		deb,fin = interval[node][0],interval[node][1]
		if time_cursor >= deb and time_cursor <= fin:
			tmp.append(node)
	tmp.sort()
	vector = "-".join(tmp)
	if vector in branch_overlap_windows:
		pass
	else:
		current_interval +=1
		branch_overlap_windows[vector] = current_interval
	time_cursor += min_period/2

# Find valid time intervals
# finds the actual time interval where ALL branches in that combination overlap simultaneously.
interval_windows={}
for branch_vector in branch_overlap_windows:
    interval_id = branch_overlap_windows[branch_vector]
    branch_list = branch_vector.split("-")
    interval_starts, interval_ends = [], []
    for branch in branch_list:
        start, end = interval[branch][0], interval[branch][1]
        interval_starts.append(start)
        interval_ends.append(end)
    interval_windows[interval_id] = [branch_list, max(interval_starts), min(interval_ends)]

# The total number of time intervals
max_intervals = max(interval_windows.keys())

########  Simulate evolution with mutations, recombination, MGEs gain and loss ######

sequence={}
sequence["root"] = concatenated_contig
for node in roots:
	parent[node] = "root"

# Data frame for substitutions
mutation_log = []
# Data frame for recombination events
recombination_log = []
# Data frame for gene gain/loss events
gain_log = []
loss_log = []

# A value to track of the sequence length
# Because the sequence length will always remain the same due to the use of '-'
# to mark the lost regions in the genome
# Thus a dict is not necessary
sequence_length = len(concatenated_contig)
# Initialize counters for total events
total_mutation,total_recomb,total_gain,total_loss = 0,0,0,0

nu_sum=[]
recomb_size_sum=[]
# Iterate through each time interval to simulate evolution
for interval_idx in range(1, max_intervals):
	interval_node = interval_windows[interval_idx][0]

	# Update sequences and the length for all nodes in the current interval
	for node_entry in interval_node:
		if node_entry in sequence:
			pass
		else:
			sequence[node_entry] = sequence[parent[node_entry]]
	
	interval_start,interval_end = interval_windows[interval_idx][1],interval_windows[interval_idx][2]
	interval_duration = interval_end - interval_start
	all_events=[]

	for node_entry in interval_node:
		# Count of mutations in total during the time segment
		# The unit of the branch length we use in this version is Year
		# The total count of mutation generate within the time period is
		# branch_length(year) * mutation rate(snp per site per year)
		# The valid sequence length is used to determine the total mutation count, instead of the full sequence length
		valid_length = sum(1 for base in sequence[node_entry] if base != "-")
		node_mutation = numpy.random.poisson(valid_length * sub_rate * interval_duration)
		for _ in range(node_mutation):
			event_entries = node_entry + "_mutation"
			all_events.append(event_entries)
		
		# Count of recombination events in total during the time segment
		node_recomb = numpy.random.poisson(recom_rate * interval_duration)
		for _ in range(node_recomb):
			# Here I introduce 2 recombination events
			# _rint is the internal recombination event, where we have donor and recipient
			# _rext is the external recombination event, where we assume there are potential undetected donor
			# the sequence from the donor will contain random snps generated using numpy.random.poisson(nu * len(recombination_sequence))
			event_entries = node_entry + numpy.random.choice(["_recomb_int", "_recomb_ext"])
			all_events.append(event_entries)

		# MGEs gains and losses
		if not use_mge or gain_rate == 0 or loss_rate == 0 or chr_bins is None or mge_fasta is None or mge_entropy is None:
			node_mge_gain, node_mge_loss = 0, 0
		else:
			node_mge_gain = numpy.random.poisson(gain_rate * interval_duration)
			for _ in range(node_mge_gain):
				all_events.append(node_entry + "_gain")

			node_mge_loss = numpy.random.poisson(loss_rate * interval_duration)
			for _ in range(node_mge_loss):
				all_events.append(node_entry + "_loss")

	# shuffle the mutation and recombination events
	random.shuffle(all_events)
	
	for event_entry in all_events:
		event_node = event_entry.split("_")[0]
		event_category = event_entry.split("_")[1]+"_"+event_entry.split("_")[2] if len(event_entry.split("_"))==3 else event_entry.split("_")[1]

		if event_category == "mutation":
			# mutated genomes and mutation info
			# Use weighted or uniform mutation explicitly based on configuration
			if use_weight_mutate:
				mutated_seq, mutation_info = mutator.weight_mutate(sequence[event_node], node=event_node)
			else:
				mutated_seq, mutation_info = mutator.uniform_mutate(sequence[event_node], node=event_node)
			# Update the sequence with the mutated version
			sequence[event_node] = mutated_seq
			# record mutation position
			mutation_log.append(mutation_info)
			total_mutation+=1

		elif event_category == "gain":
			if use_mge and gain_rate != 0 and chr_bins is not None and mge_fasta is not None and mge_entropy is not None:
				# mutated genomes and gene gain info
				mutated_seq, gene_gain_info = mgeser.mge_gain(sequence[event_node], node=event_node)
				# Update the sequence with the mutated version
				sequence[event_node] = mutated_seq
				# record gene gain event
				gain_log.append(gene_gain_info)
				total_gain += 1
			else:
				print("MGE gain simulation skipped due to incomplete data.")
		
		elif event_category == "loss":
			if use_mge and loss_rate != 0 and chr_bins is not None and mge_fasta is not None and mge_entropy is not None:
				# mutated genomes and gene loss info
				mutated_seq, gene_loss_info = mgeser.mge_loss(sequence[event_node], node=event_node)
				# Update the sequence with the mutated version
				sequence[event_node] = mutated_seq
				# Note: sequence length remains unchanged for gene loss as '-' are used to mark lost regions
				# record gene loss event
				loss_log.append(gene_loss_info)
				total_loss += 1
			else:
				print("MGE loss simulation skipped due to incomplete data.")
			
		elif event_category == "recomb_int":
			# The recombination between the same lineage happens only when there are more than 1 node(genome) in the interval
			if len(interval_node)>1:
				# the donor is randomly picked from the other genomes in the same interval
				# the recipient is the node itself
				# the recombination size is generated using exponential distribution with mean size provided
				# the start position is randomly picked from the recipient genome
				# the end position is determined by the start position and recombination size
				# if the end position exceeds the genome length, it is set to the genome length
				# the sequence from the donor is inserted into the recipient at the specified position
				# using 'nu' to generate random mutations in the recombination event
				# with this, snps can now originate from both recombination and random poisson process
				donor_candidate = [n for n in interval_node if n != event_node]
				donor = random.choice(donor_candidate)

				recomb_size = numpy.random.exponential(mean_recomb_size)
				if recomb_size < min_recomb_size:
					recomb_size = min_recomb_size
				
				start = random.choice(range(sequence_length))
				end = int(round(start + recomb_size,0))

				donor_seq = sequence[donor][start:end]
				# Check if "-" in the donor_seq
				if "-" in donor_seq:
					# Count the number of '-' characters in the recombination region.
					# Each '-' represents a gap, which will extend the recombination region by that many base pairs.
					num_gaps = donor_seq.count("-")
					end += num_gaps
				# Edge handling: Make sure end does not exceed sequence length
				if end >= sequence_length:
					end = sequence_length - 1
				# Re-gain the donor_seq after adjusting the end position
				donor_seq = sequence[donor][start:end]
				# If the old_rext_seq contains only '-', skip this recombination event
				# as it means the entire region is lost and cannot be recombined
				# This is a rare case but should be handled to avoid errors
				donor_seq_valid_len = len(donor_seq.replace("-",""))

				if donor_seq_valid_len == 0:
					continue
				else:
					receiver_seq = sequence[event_node][start:end]
					receiver_seq_valid_len = len(receiver_seq.replace("-",""))
					if receiver_seq_valid_len > 0:
						# snp_count to record the snps in the recombination events
						snp_count = 0
						snp_pos = []
						for pos_idx in range(len(receiver_seq)):
							donor_nt,receiver_nt = donor_seq[pos_idx] , receiver_seq[pos_idx]
							if donor_nt != receiver_nt:
								snp_pos.append(str(start+pos_idx))
								snp_count +=1
						nu_idx = int(snp_count) / len(donor_seq)

						total_recomb +=1
						recombination_info = {
							'donor' : str(donor),
							'recipient' : str(event_node),
							'start' : int(start),
							'end' : int(end),
							'size' : donor_seq_valid_len,
							'snp_introduced' : int(snp_count),
							'snp_position' : ";".join(snp_pos)
							}
						recombination_log.append(recombination_info)
						sequence[event_node] = sequence[event_node][:start] + sequence[donor][start:end] + sequence[event_node][end:]
						nu_sum.append(nu_idx)
						recomb_size_sum.append(donor_seq_valid_len)

		elif event_category == "recomb_ext":
			# using 'nu' to generate random mutations in the recombination event
			# with this, snps can now originate from both recombination and random poisson process
			# the first part of this is exactly the same as rint
			recomb_size = numpy.random.exponential(mean_recomb_size)
			
			if recomb_size < min_recomb_size:
				recomb_size = min_recomb_size

			start = random.choice(range(sequence_length))
			end = int(round(start + recomb_size,0))
			old_rext_seq = sequence[event_node][start:end]
			
			# Check if there is any "-" in the old_rext_seq
			if "-" in old_rext_seq:
				# Count the number of '-' characters in the recombination region.
				# Each '-' represents a gap, which will extend the recombination region by that many base pairs.
				num_gaps = old_rext_seq.count("-")
				end += num_gaps

			# Edge handling: Make sure end does not exceed sequence length
			if end >= sequence_length:
				end = sequence_length - 1

			# Re-gain the old_rext_seq after adjusting the end position
			old_rext_seq = sequence[event_node][start:end]
			# If the old_rext_seq contains only '-', skip this recombination event
			# as it means the entire region is lost and cannot be recombined
			# This is a rare case but should be handled to avoid errors
			old_rext_seq_valid_len = len(old_rext_seq.replace("-",""))
			if old_rext_seq_valid_len == 0:
				continue
			# If there are valid bases in old_rext_seq, proceed with recombination
			else:
				new_rext_seq = list(old_rext_seq)
				# Only consider positions with valid bases for SNP introduction
				rext_snp = numpy.random.poisson(nu * old_rext_seq_valid_len)
				# snp_pos to keep track of the position
				snp_pos = []
				for _ in range(rext_snp):
					# Only select positions that are not '-'
					valid_positions = [pos for pos, base in enumerate(old_rext_seq) if base != '-']
					if valid_positions:
						rext_snp_pos = numpy.random.choice(valid_positions)
						snp_pos.append(str(start + rext_snp_pos))
						rext_old_base = new_rext_seq[rext_snp_pos]
						rext_new_base = numpy.random.choice(probability[rext_old_base])
						new_rext_seq[rext_snp_pos] = rext_new_base
				# Convert back to string properly
				new_rext_seq_str = ''.join(new_rext_seq)
				# Calculate nu_idx based on the valid length of old_rext_seq
				nu_idx = int(rext_snp) / old_rext_seq_valid_len
				total_recomb += 1
				recombination_info = {
					# Now here I specify that this is an external recombination
					'donor' : "external",
					'recipient' : str(event_node),
					'start' : int(start),
					'end' : int(end),
					'size' : old_rext_seq_valid_len,
					'snp_introduced' : int(rext_snp),
					'snp_position' : ";".join(snp_pos)
					}
				
				recombination_log.append(recombination_info)
				sequence[event_node] = sequence[event_node][:start] + new_rext_seq_str + sequence[event_node][end:]
				nu_sum.append(nu_idx)
				recomb_size_sum.append(old_rext_seq_valid_len)

# Export Recombination log
recombination_log = pd.DataFrame(recombination_log)
recombination_log.to_csv(os.path.join(output_path, "recombination.csv"), index=False)

# Export Substitution log
mutation_log = pd.DataFrame(mutation_log)
mutation_log.to_csv(os.path.join(output_path, "mutations.csv"), index=False)

# Export Gene Gain/Loss log
gain_log = pd.DataFrame(gain_log)
gain_log.to_csv(os.path.join(output_path, "gain.csv"), index=False)
loss_log = pd.DataFrame(loss_log)
loss_log.to_csv(os.path.join(output_path, "loss.csv"), index=False)

# Export other logs
# Export recombination size distribution
with open(output_path + "nu.txt", "w") as h:
    h.write(f"{sum(nu_sum) / len(nu_sum)}\n")

# Read names mapping file
with open(intermediate_path + "names.txt", "r") as f:
    rename = dict(line.strip().split("\t") for line in f)

# Export the simulated genomes
# both tips and nodes
with open(output_path + "simulated_genomes.fasta", "w") as h:
	for node in sequence:
		if node != "root":
			if branch_type[node] == "tip":
				h.write(f">{rename[node]}\n{sequence[node]}\n")
			elif branch_type[node] == "branch":
				h.write(f">{node}\n{sequence[node]}\n")

# Calculate values once
rm_value = 0.0 if recom_rate == 0 else recom_rate * mean_recomb_size * sum(nu_sum) / len(nu_sum)
nu_value = "NA" if recom_rate == 0 else sum(nu_sum) / len(nu_sum)

# Export the concatenated chromosome bins
if mgeser.chr_bins is not None:
	mgeser.chr_bins_df.to_csv(os.path.join(intermediate_path, "chromosome_bins.csv"), index=False)

# Export the weight mutation DataFrame
if use_weight_mutate and mutator.weight_mutation_df is not None:
	mutator.weight_mutation_df.to_csv(os.path.join(intermediate_path, "weight_mutation.csv"), index=False)
else:
	mutator.position_df.to_csv(os.path.join(intermediate_path, "position_info.csv"), index=False)

# Export rm.txt
with open(intermediate_path + "rm.txt", "w") as h:
   h.write(f"r/m= {rm_value}\n")

# Export stats.txt
with open(output_path + "stats.txt", "w") as h:
   h.write(f"rm= {rm_value}\nNu= {nu_value}\nTotal Mutations= {total_mutation}\nTotal Recombination Events= {total_recomb}\nTotal Gain Events= {total_gain}\nTotal Loss Events= {total_loss}\n")

print(f"rm: {rm_value}\nNu: {nu_value}\nTotal Mutations: {total_mutation}\nTotal Recombination Events: {total_recomb}\nTotal Gain Events: {total_gain}\nTotal Loss Events: {total_loss}\n")

# Use snp-dist to calculate the distance matrix
print("\n" + "="*50)
print("Calculating SNP distance matrix using snp-dists")
print("="*50)
snp_fasta = f"{output_path}simulated_genomes.fasta"
snp_matrix = f"{output_path}snp_matrix.txt"
snp_dists_cmd = ["snp-dists", "-b", snp_fasta]
with open(snp_matrix, "wb") as out_log:
	proc = subprocess.run(snp_dists_cmd, stdout=out_log, stderr=subprocess.PIPE)
if proc.returncode != 0:
	print(f"Error running snp-dists (exit {proc.returncode}): {proc.stderr.decode().strip()}")
else:
	print(f"snp-dists completed, output written to {snp_matrix}")
print("SNP distance matrix calculation completed.")
