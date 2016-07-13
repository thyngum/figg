import logging
import os

"""
Core figg module
"""

import figg_dist.figg_dist as figg_dist
import figg_nj.figg_nj as figg_nj
import figg_output.figg_output as figg_output 

def run_figg(input_file, is_circular, output_format, verbose):

	# Global names
	genome_labels = []                # List of genome labels
	gene_orders = []                  # List with the order of genes in each input genome
	num_genomes = 0                   # Number of genomes in input
	ref_order = []		              # Reference gene order
	ref_matrix = []                   # Reference adjacency matrix
	first_seen = []                   # Genome in which each gene was first seen
	''' >>> TODO <<<: Convert first_seen to a dictionary instead of a positional list? '''
	adj_matrices = []                 # Set of adjacency matrices for each genome in input
	dist_matrix = []                  # Matrix of observed differences
	corrected_dist_matrix = []        # Matrix of corrected differences

	# Read input file
	f = open(input_file, "rU")
	temp = f.read()
	temp = temp.split('>')
	del temp[0]

	# Parse the input
	for i in range(len(temp)):
		genome_labels.append(temp[i].split('\n')[0])
		gene_orders.append(temp[i].split('\n')[1])
	for i in range(len(gene_orders)):
		gene_orders[i] = gene_orders[i].split(" ")
	num_genomes = len(gene_orders)

	# If genomes are circular then append the first gene at the end of each gene list
	if is_circular:
		[i.append(i[0]) for i in gene_orders]

	# Initialize the reference order with that of the first genome and delete it from the list 
	''' >>> TODO <<<: Perhaps find a way to avoid deleting it from the list '''
	ref_order = gene_orders[0] 
	del gene_orders[0]
	first_seen = [genome_labels[0]]*len(ref_order)

	# If genomes are circular the first gene will be duplicated at the end of first_seen
	''' >>> TODO <<<: This issue will be solved if first_seen is converted to a dict  '''
	if ( is_circular ):
		del first_seen[-1]
	
	# Initialize the reference matrix
	num_genes = len(set(ref_order))
	index_list = [ref_order.index(i) for i in ref_order]
	ref_matrix = [[0]*num_genes for i in range(num_genes)]
	for i in range(len(index_list) - 1):
		if "-" in ref_order[i + 1]:
			ref_matrix[index_list[i]][index_list[i + 1]] = -1
		else:
			ref_matrix[index_list[i]][index_list[i + 1]] = 1

	# Extend the reference order and matrix by looking for new genes in all other input genomes
	for i in range(1,num_genomes):
		for j in gene_orders[i - 1]: 
			a = j not in ref_order
			b = j.replace('-','') not in ref_order
			if a & b:
				ref_order.append(j)
				ref_matrix.append([0]*len(ref_matrix))
				first_seen.append(genome_labels[i])
				for x in range(len(ref_matrix)):
					ref_matrix[x] += [0]
	ref_order.pop(ref_order[1:].index(ref_order[0]) + 1)

	# Compute all the matrices
	adj_matrices = figg_dist.adj_matrix_set(ref_matrix, ref_order, gene_orders)
	dist_matrix = figg_dist.dist_matrix(adj_matrices)
	pos_freq_matrix = figg_dist.freq_matrix_pos(adj_matrices)
	neg_freq_matrix = figg_dist.freq_matrix_neg(adj_matrices)
	corrected_dist_matrix = figg_dist.dist_matrix_corrected(adj_matrices, pos_freq_matrix, neg_freq_matrix)

	# Build the NJ tree 
	tree, heights = figg_nj.nj(dist_matrix, genome_labels, [])
	tree_corrected, heights_corrected = figg_nj.nj(corrected_dist_matrix, genome_labels, [])

	# Write the output
	abs_path = os.path.abspath(input_file)
	output_path, input_filename = os.path.split(abs_path)
	output_prefix = output_path + '/' + '.'.join(input_filename.split('.')[:-1])

	if ( output_format == 'text' or output_format == 'all' ):
		figg_output.print_matrix_to_file(dist_matrix, output_prefix + "_dist_matrix.tsv", genome_labels)	
		figg_output.print_matrix_to_file(corrected_dist_matrix, output_prefix + "_corrected_dist_matrix.tsv", genome_labels)

	if ( output_format == 'mega' or output_format == 'all' ):
		figg_output.print_mega_format(dist_matrix, genome_labels, output_prefix + "_dist_matrix.meg")	
		figg_output.print_mega_format(corrected_dist_matrix, genome_labels, output_prefix + "_corrected_dist_matrix.meg")

	if ( output_format == 'phylip' or output_format == 'all' ):
		figg_output.print_phylip_format(dist_matrix, genome_labels, output_prefix + "_dist_matrix.phy")	
		figg_output.print_phylip_format(corrected_dist_matrix, genome_labels, output_prefix + "_corrected_dist_matrix.phy")	

	if ( output_format == 'nexus' or output_format == 'all' ):
		figg_output.print_nexus_format(dist_matrix, genome_labels, output_prefix + "_dist_matrix.nex")	
		figg_output.print_nexus_format(corrected_dist_matrix, genome_labels, output_prefix + "_corrected_dist_matrix.nex")

	if ( verbose ):

		# Print basic info
		print "Number of genomes: %i" % num_genomes
		print "Number of genes in the workspace: %i" % len(ref_order)
		print "  Gene\tFirst seen in"
		for i in range(len(ref_order)):
			print "  %s\t%s" % (ref_order[i], first_seen[i])
		print ""

		# Print adjacency matrices for each genome
		print "Adjacency matrices for each genome:\n"
		for i in range(num_genomes):
			figg_output.print_matrix(adj_matrices[i], ref_order)		

		# Print frequency matrices
		print "Positive and negative frequency matrices:\n"
		figg_output.print_matrix(pos_freq_matrix, ref_order)
		figg_output.print_matrix(neg_freq_matrix, ref_order)

		# Print distance matrices
		print "Uncorrected distance matrix:\n"		
		figg_output.print_matrix(dist_matrix, genome_labels)
		print "Corrected distance matrix:\n"
		figg_output.print_matrix(corrected_dist_matrix, genome_labels)

		# Print the NJ tree
		print "NJ tree from uncorrected distance matrix:\n"
		print tree, "\n"
		print "NJ tree from corrected distance matrix:\n"
		print tree_corrected, "\n"

		print "Output saved to '%s'." % output_path
