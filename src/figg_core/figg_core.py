import logging

import figg_dist.figg_dist as figg_dist
import figg_nj.figg_nj as figg_nj
import figg_output.figg_output as figg_output 

def run_figg(input_file, is_circular, output_format, verbose):

	# Global names
	genome_labels = []              # List of genome labels
	gene_orders = []                # List with the order of genes in each input genome
	num_genomes = 0                 # Number of genomes in input
	reference_order = []			# Reference gene order
	reference_matrix = []           # Reference adjacency matrix
	first_seen = []                 # Genome in which each gene was first seen

	# Reads input file
	f = open(input_file, "rU")
	temp = f.read()
	temp = temp.split('>')
	del temp[0]

	# Parses the input
	for i in range(len(temp)):
		genome_labels.append(temp[i].split('\n')[0])
		gene_orders.append(temp[i].split('\n')[1])
	for i in range(len(gene_orders)):
		gene_orders[i] = gene_orders[i].split(" ")
	num_genomes = len(gene_orders)

	# If genomes are circular then appends the first gene at the end of each gene order
	if is_circular:
		[i.append(i[0]) for i in gene_orders]

	# Initializes the reference order with that of the first genome and delete it from the list 
	reference_order = gene_orders[0] 
	del gene_orders[0]
	first_seen = [genome_labels[0]]*len(reference_order)
	if ( is_circular ):
		del first_seen[-1]
	
	# Initializes the reference matrix
	ref_index_list = [reference_order.index(i) for i in reference_order]
	reference_matrix = [[0]*len(set(reference_order)) for i in range(len(set(reference_order)))]
	for i in range(len(ref_index_list) - 1):
		if "-" in reference_order[i + 1]:
			reference_matrix[ref_index_list[i]][ref_index_list[i + 1]] = -1
		else:
			reference_matrix[ref_index_list[i]][ref_index_list[i + 1]] = 1

	# Extends the reference order and matrix by looking for new genes in all other input genomes
	k = 1
	for i in gene_orders:
		for j in i:
			a = j not in reference_order
			b = j.replace('-','') not in reference_order
			if a & b:
				reference_order.append(j)
				reference_matrix.append([0]*len(reference_matrix))
				first_seen.append(genome_labels[k])
				for x in range(len(reference_matrix)):
					reference_matrix[x] += [0]
       		k += 1
	rem = reference_order[1:].index(reference_order[0]) + 1
	reference_order.pop(rem)

	if ( verbose ):
		print "Number of genomes: %i" % num_genomes
		print "Number of genes in the workspace: %i" % len(reference_order)
		print "  Gene\tFirst seen in"
		for i in range(len(reference_order)):
			print "  %s\t%s" % (reference_order[i], first_seen[i])
		print ""

	distance_matrix, M = figg_dist.dmatrix(reference_matrix, reference_order, gene_orders)
	corrected_distance_matrix, positive_freq_matrix, negative_freq_matrix = figg_dist.cdmatrix(M)

	if ( verbose ):

		# Print adjacency matrices for each genome
		print "Adjacency matrices for each genome:\n"
		for i in range(num_genomes):
			figg_output.print_matrix(M[i], reference_order)		

		# Print frequency matrices
		print "Positive and negative frequency matrices:\n"
		figg_output.print_matrix(positive_freq_matrix, reference_order)
		figg_output.print_matrix(negative_freq_matrix, reference_order)

		# Print distance matrices
		print "Uncorrected distance matrix:\n"		
		figg_output.print_matrix(distance_matrix, genome_labels)
		print "Corrected distance matrix:\n"
		figg_output.print_matrix(corrected_distance_matrix, genome_labels)

	# Write the output
	# figg_output.print_matrix_to_file(distance_matrix, "distance_matrix.tsv", genome_labels)	# For these, find how to get the path of input files
	# figg_output.print_matrix_to_file(corrected_distance_matrix, "corrected_distance_matrix.tsv", genome_labels)	

	# figg_output.print_matrix_to_file(positive_freq_matrix, reference_order, "fplus_admat.txt")
	# figg_output.print_matrix_to_file(negative_freq_matrix, reference_order, "fneg_admat.txt")
	# figg_output.print_mega_format(corrected_distance_matrix, genome_labels, "corrected_distmat.meg")
	
	# Build the NJ tree 
	# n_CD = [corrected_distance_matrix[i][:] for i in range(len(corrected_distance_matrix))]
	# nj(n_CD,refs,[])
