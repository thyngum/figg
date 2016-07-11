import logging

"""
Functions to work with distance matrices
"""

def AD(matrix_1, matrix_2):
	"Computes the measure of differences between two adjacency matrices, AD"

	d = 0 
	for i in range(len(matrix_1)):
		for j in range(len(matrix_2)):
			if matrix_1[i][j] == matrix_2[i][j]:
				d += 1
	
	return(d)


def workspace(rseq, aseq):
	n = len(set(rseq))
	K = [[0]*n for i in range(n)]

	ref = [rseq.index(i) if i in rseq else rseq.index(i.replace('-','')) for i in aseq]
	for i in range(len(ref)-1):
		if "-" in list(aseq[i+1]):
			K[ref[i]][ref[i+1]] = -1
		else:
			K[ref[i]][ref[i+1]] = 1 
	return K


def dmatrix(ref_string_adjacency_matrix, ref_string, strings):
	"Computes the uncorrected distance matrix (matrix of AD)" 

	num_genomes = len(strings) + 1

	temp_matrix = []
	temp_matrix.append(ref_string_adjacency_matrix)
	[temp_matrix.append(workspace(ref_string, strings[i])) for i in range(len(strings))]
	distance_matrix = [[0]*num_genomes for i in range(num_genomes)]
	for i in range(num_genomes):
		for j in range(i):
			distance_matrix[i][j] = AD(temp_matrix[i], temp_matrix[j])
			
	return distance_matrix, temp_matrix 


def fpos(temp_matrix):
	"Computes the positive frequency matrix"

	n = len(temp_matrix[0])
	K = [[0]*n for i in range(n)]
	for i in range(len(K)):
		for j in range(len(K)):
			K[i][j] = round(sum([1 for x in range(len(temp_matrix)) if temp_matrix[x][i][j] == 1]) / float(len(temp_matrix)),4)
	
	return K


def fneg(temp_matrix):
	"Computes the negative frequency matrix"

	n = len(temp_matrix[0])
	K = [[0]*n for i in range(n)]
	for i in range(len(K)):
		for j in range(len(K)):
			K[i][j] = round(sum([1 for x in range(len(temp_matrix)) if temp_matrix[x][i][j] == -1]) / float(len(temp_matrix)),4)
	
	return K


def cdiff(i, j, matrix_1, matrix_2, positive_freqs, negative_freqs):
	"Similar to diff() but takes into account positive and negative frequencies"

	if matrix_1[i][j] == matrix_2[i][j]:
		return(0)
	elif (matrix_1[i][j] == matrix_2[i][j]) & (matrix_1[i][j] > 0 & matrix_2[i][j] > 0):
		return(1 - positive_freqs[i][j])
	elif (not(matrix_1[i][j] == matrix_2[i][j]) & (matrix_1[i][j]) < 0 | matrix_2[i][j] < 0):
		return(1 - negative_freqs[i][j])


def ADc(matrix_1, matrix_2, positive_freqs, negative_freqs):
	"Computes the corrected measure of differences between two adjacency matrices, ADc"

	dc = 0 
	for i in range(len(matrix_1)):
		for j in range(len(matrix_1)):
			dc += cdiff(i, j, matrix_1, matrix_2, positive_freqs, negative_freqs)
	
	return(dc)


def cdmatrix(temp_matrix):
	"Computes the corrected distance matrix (matrix of AD)"

	positive_freqs = fpos(temp_matrix)
	negative_freqs = fneg(temp_matrix)
	k = len(temp_matrix)
	CD =  [[0]*k for i in range(k)]
	for i in range(k):
		for j in range(i):
			CD[i][j] = ADc(temp_matrix[i],temp_matrix[j], positive_freqs, negative_freqs)
			
	return CD, positive_freqs, negative_freqs   
