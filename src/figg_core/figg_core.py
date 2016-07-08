import logging

import figg_dist.figg_dist as figg_dist
import figg_nj.figg_nj as figg_nj

def run_figg(input_file, is_circular, output_format):

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

	logging.info("Number of genomes: %i" % num_genomes)
	logging.info("Workspace dimension: %ix%i\n" % (workspace_dimension[0], workspace_dimension[0]))

	distance_matrix, M = figg_dist.dmatrix(ref_string_adj_matrix, ref_string, strings)
	corrected_distance_matrix, positive_freq_matrix, negative_freq_matrix = figg_dist.cdmatrix(M)
	
	# Print adjacency matrices for each genome
	# for i in range(num_genomes):
	# 	print_matrix(M[i], ref_string)

	# Print frequency matrices
	# print_matrix(positive_freq_matrix, ref_string)
	# print_matrix(negative_freq_matrix, ref_string)

	# Print distance matrices
	print_matrix(distance_matrix, names)
	print_matrix(corrected_distance_matrix, names)
	print_matrix_to_file(distance_matrix, "distance_matrix.tsv", names)	# For these, find how to get the path of input files
	print_matrix_to_file(corrected_distance_matrix, "corrected_distance_matrix.tsv", names)	

	"""
	print_matrix_to_file(positive_freq_matrix, ref_string, "fplus_admat.txt")
	print_matrix_to_file(negative_freq_matrix, ref_string, "fneg_admat.txt")
	print_mega_format(corrected_distance_matrix, names, "corrected_distmat.meg")
	"""
	
	# Build the NJ tree 
	"""
	n_CD = [corrected_distance_matrix[i][:] for i in range(len(corrected_distance_matrix))]
	nj(n_CD,refs,[])
	"""



def print_matrix(matrix, labels = False):
	"Prints a matrix with optional row labels"
	
	string = ""
	k = 0
	for i in matrix:
		n = 0
		for j in i:
			if n == 0:
				if ( labels ):
					string += labels[k] + '\t' + '%.4f'%j
				else:
					string += '%.4f'%j
			else:
				string += '\t' + '%.4f'%j
			n += 1
		string += '\n'
		k += 1

	print string


def print_matrix_to_file(matrix, filename, labels = False):

	f = open(filename,'w')
	string = ""
	k = 0
	for i in matrix:
		n = 0
		for j in i:
			if n == 0:
				if ( labels ):
					string += labels[k] + '\t' + '%.4f'%j
				else:
					string += '%.4f'%j
			else:
				string += '\t' + '%.4f'%j
			n += 1
		string += '\n'
		k += 1
	f.write(string)
	f.close()

""""
def print_mega_format(matrix, labels, filename):

	num_genomes = len(matrix)

	f = open(filename, 'w')

	max_num_length = 0
	for i in matrix:
		if (len(str(int(max(i)))) > max_num_length):
			max_num_length = len(str(int(max(i))))

	# Header
	text = ("#mega\n!Title: Corrected adjacency distance matrix;\n" +
		"!Format DataType=Distance DataFormat=LowerLeft NTaxa=" + str(num_genomes) + ";\n\n")

	# List of labels
	for i in range(num_genomes):
		text += '[' + " "*(len(str(num_genomes)) - len(str(i + 1))) + str(i + 1) + '] ' + "#" + labels[i] +"\n"
	
	# First row of distance matrix
	text += '\n[' + '\t'
	x = [str(i) for i in range(1, num_genomes + 1)]
	for q in x[:-1]: 
		text += q + '\t'
	text += x[-1] + ']\n'

	# Body of distance matrix
	k = 1
	text += '[ 1]'
	for i in range(num_genomes):
		n = 0
		for j in range(i):
			if n == 0:
				text += '[' + " "*(len(str(num_genomes)) - len(str(i + 1))) + str(i + 1) + ']\t' + '%.4f'%(matrix[i][j])
			else:
				#g = max_num_length - len(str(int(matrix[i][j])))
				text += '\t' + '%.4f'%(matrix[i][j])
			n += 1
		text = text + '\n'
		k += 1

	f.write(text)
	f.close()
"""
