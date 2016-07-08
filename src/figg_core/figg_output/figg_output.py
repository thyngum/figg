import logging

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
	"Prints a matrix to a tab-separated (TSV) file with optional row labels"

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
