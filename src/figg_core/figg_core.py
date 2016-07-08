import logging

import figg_dist.figg_dist as figg_dist
import figg_nj.figg_nj as figg_nj
import figg_output.figg_output as figg_output 

def run_figg(input_file, is_circular, output_format, verbose):

	# Global names
	names = []                      # List of genome labels
	strings = []                    # List of gene strings
	num_genomes = 0                 # Number of genomes in input
	ref_string = []                 # Reference gene string   
	ref_string_adj_matrix = []      # Adjacency matrix for the reference string (a.k.a. the workspace)
	workspace_dimension = 0         # Number of genes in the workspace

	# Reads input file
	f = open(input_file, "rU")
	temp = f.read()
	temp = temp.split('>')
	del temp[0]

	# Get labels and gene strings
	for i in range(len(temp)):
		names.append(temp[i].split('\n')[0])
		strings.append(temp[i].split('\n')[1])
	for i in range(len(strings)):
		strings[i] = strings[i].split(" ")
	num_genomes = len(strings)

	# If genome is circular then appends the first gene at the end of the gene list
	if is_circular:
		[i.append(i[0]) for i in strings]

	# Sets the reference string and delete it from the list of strings 
	ref_string = strings[0] 
	del strings[0] 
	first_seen = [names[0]]*len(ref_string)

	# Computes the adjacency matrix for the reference string
	n = len(set(ref_string))
	ref = [ref_string.index(i) for i in ref_string]
	ref_string_adj_matrix = [[0]*n for i in range(n)]
	for i in range(len(ref) - 1):
		if "-" in ref_string[i + 1]:
			ref_string_adj_matrix[ref[i]][ref[i + 1]] = -1
		else:
			ref_string_adj_matrix[ref[i]][ref[i + 1]] = 1

	# Goes through all other sequences and looks for new genes. Extends the dimension 
	# of the workspace by adding columns and rows. After this, reference string will include  
	# all found genes, with each new gene added to its end. 
	k = 1
	for i in strings:
		for j in i:
			a = j not in ref_string
			b =  j.replace('-','') not in ref_string
			if a & b:
				ref_string.append(j)
				ref_string_adj_matrix.append([0]*len(ref_string_adj_matrix))
				first_seen.append(names[k])
				for x in range(len(ref_string_adj_matrix)):
					ref_string_adj_matrix[x] += [0]
		k += 1
	rem = ref_string[1:].index(ref_string[0]) + 1
	ref_string.pop(rem)

	workspace_dimension = [len(ref_string_adj_matrix), len(ref_string_adj_matrix[0])]

	if ( verbose ):
		print "Number of genomes: %i" % num_genomes
		print "Workspace dimension: %ix%i\n" % (workspace_dimension[0], workspace_dimension[0])

	distance_matrix, M = figg_dist.dmatrix(ref_string_adj_matrix, ref_string, strings)
	corrected_distance_matrix, positive_freq_matrix, negative_freq_matrix = figg_dist.cdmatrix(M)
	
	if ( verbose ):

		# Print adjacency matrices for each genome
		print "Adjacency matrices for each genome:\n"
		for i in range(num_genomes):
			figg_output.print_matrix(M[i], ref_string)		

		# Print frequency matrices
		print "Positive and negative frequency matrices:\n"
		figg_output.print_matrix(positive_freq_matrix, ref_string)
		figg_output.print_matrix(negative_freq_matrix, ref_string)

		# Print distance matrices
		print "Uncorrected distance matrix:\n"		
		figg_output.print_matrix(distance_matrix, names)
		print "Corrected distance matrix:\n"
		figg_output.print_matrix(corrected_distance_matrix, names)

	# Write the output
	figg_output.print_matrix_to_file(distance_matrix, "distance_matrix.tsv", names)	# For these, find how to get the path of input files
	figg_output.print_matrix_to_file(corrected_distance_matrix, "corrected_distance_matrix.tsv", names)	

	"""
	figg_output.print_matrix_to_file(positive_freq_matrix, ref_string, "fplus_admat.txt")
	figg_output.print_matrix_to_file(negative_freq_matrix, ref_string, "fneg_admat.txt")
	figg_output.print_mega_format(corrected_distance_matrix, names, "corrected_distmat.meg")
	"""
	
	# Build the NJ tree 
	"""
	n_CD = [corrected_distance_matrix[i][:] for i in range(len(corrected_distance_matrix))]
	nj(n_CD,refs,[])
	"""
