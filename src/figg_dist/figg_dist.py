import logging

"""
Functions to work with distance matrices
"""

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


def bij(i, j, m, n):
	# m and n are matrices
	# i is the row position and j is the column position
	if m[i][j] == n[i][j]:
		return(0)
	else:
		return(1)


def AD(m, n):
	d = 0 
	for i in range(len(m)):
		for j in range(len(m)):
			d = d + bij(i,j,m,n)
	return(d)


def dmat(R, rseq, seqs):
	# rseq is the reference string
	# seqs is a lis of the other sequences
	k = len(seqs)+1
	N = []
	n = len(set(rseq))
	# R is the adjacency matrix of the reference sequence for which we use to 
	# develop the workspace 
	M = []
	M.append(R)
	[M.append(workspace(rseq,seqs[x])) for x in range(len(seqs))]
	D =  [[0]*k for i in range(k)]
	for i in range(k):
		for j in range(i):
			D[i][j] = AD(M[i],M[j])
			
	return D,M 


def fplus(M):
	n = len(M[0])
	K = [[0]*n for i in range(n)]
	for i in range(len(K)):
		for j in range(len(K)):
			K[i][j] = round(sum([1 for x in range(len(M)) if M[x][i][j] == 1])/float(len(M)),4)
	
	return K


def fneg(M):
	n = len(M[0])
	K = [[0]*n for i in range(n)]
	for i in range(len(K)):
		for j in range(len(K)):
			K[i][j] = round(sum([1 for x in range(len(M)) if M[x][i][j] == -1])/float(len(M)),4)
	
	return K


def cij(i, j, m, n, fp, fn):
	# m and n are matrices
	# i is the row position and j is the column position
	if m[i][j] == n[i][j]:
		return(0)
	elif (m[i][j] == n[i][j])&(m[i][j] > 0 & n[i][j] > 0):
		return(1-fp[i][j])
	elif (not(m[i][j] == n[i][j])&(m[i][j]) < 0 | n[i][j] < 0):
		return(1-fn[i][j])


def ADc(m, n, fp, fn):
	dc = 0 
	for i in range(len(m)):
		for j in range(len(m)):
			dc = dc + cij(i,j,m,n,fp,fn)
	return(dc)

def Cdmat(M):
	fp = fplus(M)
	fn = fneg(M)
	k = len(M)
	CD =  [[0]*k for i in range(k)]
	for i in range(k):
		for j in range(i):
			CD[i][j] = ADc(M[i],M[j],fp,fn)
			
	return CD,fp,fn   
