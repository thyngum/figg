import logging

import figg_dist.figg_dist as figg_dist
import figg_nj.figg_nj as figg_nj

def run_figg(input_file, is_circular, output_format):

	# Global names
	names = []						# List of genome labels
	strings = []					# List of gene strings
	num_genomes = 0					# Number of genomes in input
	ref_string_adj_matrix = []		# Adjacency matrix for the reference string (a.k.a. the workspace)

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

	# If genome is circular then appends the first gene at the end of gene list
	if is_circular:
		[i.append(i[0]) for i in strings]

	# Sets the reference string and delete it from the list of strings 
	rseq = strings[0] # rseq is the reference string
	del strings[0] # seqs is a list of the other sequences
	seqs = strings
	first_seen = [names[0]]*len(rseq)
	Rnames = tuple(rseq) # Labels in the reference string
	# k = len(seqs) + 1

	# Computes the adjacency matrix for the reference string
	n = len(set(rseq))
	ref = [rseq.index(i) for i in rseq]
	ref_string_adj_matrix = [[0]*n for i in range(n)]
	for i in range(len(ref)-1):
		if "-" in rseq[i+1]:
			ref_string_adj_matrix[ref[i]][ref[i+1]] = -1
		else:
			ref_string_adj_matrix[ref[i]][ref[i+1]] = 1

	# Goes through all other sequences in search for new genes and extends the dimension 
	# of the workspace by adding columns and rows. After this, reference will include all 
	# genes, with each new gene added to the its end. 
	k = 1
	for i in strings:
		for j in i:
			a = j not in rseq
			b =  j.replace('-','') not in rseq
			if a & b:
				rseq.append(j)
				ref_string_adj_matrix.append([0]*len(ref_string_adj_matrix))
				first_seen.append(names[k])
				for x in range(len(ref_string_adj_matrix)):
					ref_string_adj_matrix[x] += [0]
		k += 1

	rem = rseq[1:].index(rseq[0]) + 1
	rseq.pop(rem)

	dim = [len(ref_string_adj_matrix), len(ref_string_adj_matrix[0])]
	s = len(seqs) + 1	
	m = len(str(len(ref_string_adj_matrix)))

	logging.info("Number of genomes: %i" % num_genomes)
	logging.info("Workspace dimension: %ix%i\n" % (dim[0], dim[0]))

	D,M = figg_dist.dmat(ref_string_adj_matrix, rseq, seqs)
	CD,fp,fn = figg_dist.Cdmat(M)
	
	names = tuple(names)
	refs = list(names)

	n_CD = [CD[i][:] for i in range(len(CD))]
	#nj(n_CD,refs,[])

	# Print adjacency matrices for each genome
	"""
	for i in range(num_genomes):
		print_matrix(M[i], rseq)
	"""

	# Print frequency matrices
	"""
	print_matrix(fp, rseq)
	print_matrix(fn, rseq)
	"""

	# Print distance matrices
	print_matrix(D, names)
	print_matrix(CD, names)

"""
	print_matrix_to_file(fp, rseq, "fplus_admat.txt")
	print_matrix_to_file(fn, rseq, "fneg_admat.txt")

	print_matrix_to_file(CD, names, "test")
	print_mega_format(CD, names, "corrected_distmat.meg")
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


"""
def print_matrix_to_file(matrix, labels, filename):

	f = open(filename,'w')
	text = ""
	k = 0
	for i in matrix:
		n = 0
		for j in i:
			if n == 0:
				text += labels[k] + '\t' + '%.4f'%j
			else:
				text += '\t' + '%.4f'%j
			n += 1
		text += '\n'
		k += 1
	f.write(text)
	f.close()

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
