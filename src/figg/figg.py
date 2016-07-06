import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    print " Error: Please specify file!"
    quit()


infile = sys.argv[1]

say_yes = ['y', 'Y', 'yes', 'YES', 'Yes']

# asks whether the genome is circular or not
is_circular = raw_input("Are the genomes circular? [y/n]: ") in say_yes

# Reads fasta and stores in memory
try:
    fasta_sequences = SeqIO.parse(open(infile),'fasta')
    names = []
    strings = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        names.append(name)
        string = sequence.split(" ")
        if is_circular:
            string.append(string[0])
        strings.append(string)
except ValueError, e:
    print "Input format error: '%s'" % str(e)



rseq = strings[0]
del strings[0]
seqs = strings
first_seen = [names[0]]*len(rseq)

# rseq is the reference string
# seqs is a list of the other sequences
Rnames = tuple(rseq)
k = len(seqs) + 1
	
n = len(set(rseq))
# R is the adjacency matrix of the reference sequence for which we use to 
# develop the workspace 
ref = [rseq.index(i) for i in rseq]
R = [[0]*n for i in range(n)]
for i in range(len(ref)-1):
	if "-" in rseq[i+1]:
		R[ref[i]][ref[i+1]] = -1
	else:
		R[ref[i]][ref[i+1]] = 1
# Now we search through all other sequences for new genes and extend the dimension 
# of the workspace by adding columns and rows.
# Our reference sequence is now extended to include every gene, with each new each added
#to the tail of the reference string. 
k = 1
for i in strings:
	for j in i:
		a = j not in rseq
		b =  j.replace('-','') not in rseq
		if a&b:
			rseq.append(j)
			R.append([0]*len(R))
			first_seen.append(names[k])
			for x in range(len(R)):
				R[x] = R[x] + [0]
	k += 1

rem =  rseq[1:].index(rseq[0]) + 1
rseq.pop(rem)

# ------

def workspace(rseq,aseq):
	n = len(set(rseq))
	K = [[0]*n for i in range(n)]

	ref = [rseq.index(i) if i in rseq else rseq.index(i.replace('-','')) for i in aseq]
	for i in range(len(ref)-1):
		if "-" in list(aseq[i+1]):
			K[ref[i]][ref[i+1]] = -1
		else:
			K[ref[i]][ref[i+1]] = 1 
	return K




def bij(i,j,m,n):
	# m and n are matrices
	# i is the row position and j is the column position
	if m[i][j] == n[i][j]:
		return(0)
	else:
		return(1)
	



def AD(m,n):
	d = 0 
	for i in range(len(m)):
		for j in range(len(m)):
			d = d + bij(i,j,m,n)
	return(d)




def dmat(rseq,seqs):
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

def cij(i,j,m,n,fp,fn):
	# m and n are matrices
	# i is the row position and j is the column position
	if m[i][j] == n[i][j]:
		return(0)
	elif (m[i][j] == n[i][j])&(m[i][j] > 0 & n[i][j] > 0):
		return(1-fp[i][j])
	elif (not(m[i][j] == n[i][j])&(m[i][j]) < 0 | n[i][j] < 0):
		return(1-fn[i][j])




def ADc(m,n,fp,fn):
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

# ----


dim = [len(R), len(R[0])]
s = len(seqs) + 1
print "Genomes: " + str(s)
print '^       ' + 'Label'
m = len(str(len(R)))
for i in range(len(names)):
	print '^(' + ' '*(m-len(str(i+1))) + str(i+1) + ') ' + '\t' + names[i]

print "Workspace dimension: " + str(dim[0]) + " X " + str(dim[0])
print '^     ' + 'Item' + '\t\t' + 'first found in'
m = len(str(s))
for i in range(len(rseq)):
	h = len(rseq[i])
	g = max([len(k) for k in rseq])
	print '^(' + ' '*(m - len(str(i+1))) + str(i+1) + ')' + ' '*(g-len(rseq[i]) +2) + \
	rseq[i] +  ' '*(len(str(i))-1) + '\t\t\t' + '%s' % (first_seen[i])




y = 'y'
x = raw_input('Proceed with this dataset? [y/n]: ')


if x == 'y':
	D,M = dmat(rseq,seqs)
	CD,fp,fn = Cdmat(M)
	# Neighborjoining algorithm
	def nj(d,refs,heights):
		if len(d) == 2:
			heights += [d[1][0]/2]
			edge_i = d[1][0] 
 			edge_j = d[1][0] 
 			print "The neighbor adjoining algorith produced the following tree: \n"
			clust =  "(" + refs[0] + ": " + str(edge_j) + "," + refs[1] + ": " + str(edge_i) + ")"
			
			print clust + "\n"
			print "With the following node heights: \n"
			print heights , "\n\n"
			filename = sys.argv[1].split(".")[0]
			f = open(filename + ".nwk",'w')
			f.write(clust)
			f.close()

		else:
			L = len(d) - 2
			r = [sum([d[i][j] if j < i else d[j][i] for j in range(len(d))])/L for i in range(len(d))]

			MIN = d[1][0] - (r[1] + r[0])
 			index = (1,0) 
 			for i in range(2,len(d)):
				temp = min([(d[i][j]  - (r[i] + r[j])) for j in range(i)])
 				if MIN > temp:
 					MIN = temp
 					index = (i,[(d[i][j] - (r[i] + r[j])) for j in range(i)].index(MIN))
 			edge_i = .5*(d[index[0]][index[1]] + r[index[0]] - r[index[1]])
 			edge_j = .5*(d[index[0]][index[1]] + r[index[1]] - r[index[0]])

 			clust = "(" + str(refs[index[1]]) + ": " + str(edge_j) +\
 			 ", " + str(refs[index[0]]) + ": " + str(edge_i) + ")"		
  			
  			refs.remove(refs[index[0]])
  			refs.remove(refs[index[1]])
  			refs.insert(0,clust)
  			i0 = index[0]
  			i1 = index[1]
  			heights += [d[index[0]][index[1]]/2]
  			temp_d =  [(.5*(d[i][i1]+ d[i][i0])) if (i0 <i and i1 <i) else \
  			(.5*(d[i1][i]+ d[i][i0])) if (i0 <i and i1 > i) else \
  			(.5*(d[i][i1]+ d[i0][i])) if (i0 >i and i1 <i) else \
  			(.5*(d[i1][i]+ d[i0][i])) for i in range(len(d)) if i  not in index]
  			
  			d.remove(d[index[0]])
			d.remove(d[index[1]])
			for i in range(len(d)):
				d[i].remove(d[i][index[0]])
				d[i].remove(d[i][index[1]])

			for i in range(len(d)):
				d[i] = [temp_d[i]] + d[i]

			d.insert(0,[0 for i in range(len(d[0]))])
			return(nj(d,refs,heights))
	
	names = tuple(names)
	refs = list(names)
	
	n_CD = [CD[i][:] for i in range(len(CD))]
	nj(n_CD,refs,[])
	# writing matrices to files
	ext =  raw_input('Select output format: \n[M]EGA (.meg)\n[P]hylip (.phy)\n[N]exus (.nex)\n[A]ll\n[No]ne\n')

	admats = [0]*len(M)
	for i in range(len(M)):
		admats[i] = names[i] + '_admat.txt' 


	if  ((ext == 'M')|(ext == 'A')):
		distmat = 'uncorrected_distmat.meg'
		corr_dm = 'corrected_distmat.meg'
		f = open(distmat,'w')
		t = ("#mega\n!Title: Uncorrected adjacency distance matrix;\n" +
			"!Format DataType=Distance DataFormat=LowerLeft NTaxa=" + str(s) + ";\n\n")
		for i in range(s):
			if len(str(i+1)) == 1:
				t = t + '[ ' + str(i+1) + '] ' + "#" + names[i] +"\n"
			else:
				t = t +'[' + str(i+1) + '] ' + "#" + names[i] +"\n"
		t += '\n[' + '\t\t'
		x = [str(i) for i in range(1,s+1)]
		for q in x[:-1]: t+=q + '\t'
		t += x[-1] + ']\n'
		
		k = 1
		t += '[ ' + str(1) + ']' + '\t  ' 
		for i in range(len(D)):
			n = 0
			for j in range(i):
				if n == 0:
					if len(str(k)) == 1:
						t += '[ ' + str(k) + ']' + '\t\t' + '%s'%(D[i][j])
					else:
						t += '[' + str(k) + ']' + '\t\t' + '%s'%(D[i][j])
				else:
					t = t + '\t' + '%s'%(D[i][j])
				n = n +1
			t = t + '\n'
			k +=1
		f.write(t)
		f.close()

		f = open(corr_dm,'w')
		MAX = 0
		for i in CD:
			if (len(str(int(max(i)))) > MAX):
				MAX = len(str(int(max(i))))
		t = ("#mega\n!Title: Corrected adjacency distance matrix;\n" +
			"!Format DataType=Distance DataFormat=LowerLeft NTaxa=" + str(s) + ";\n\n")
		for i in range(s):
			if len(str(i+1)) == 1:
				t = t + '[ ' + str(i+1) + '] ' + "#" + names[i] +"\n"
			else:
				t = t +'[' + str(i+1) + '] ' + "#" + names[i] +"\n"
		t += '\n[' + '\t\t'
		x = [str(i) for i in range(1,s+1)]
		for q in x[:-1]: t+=q + '\t    '
		t += x[-1] + ']\n'
		k = 1
		t += '[ ' + str(1) + ']' + '\t    ' 
		for i in range(len(CD)):
			n = 0
			for j in range(i):
				if n == 0:
					if len(str(k)) == 1:
						t += '[ ' + str(k) + ']' + '\t      ' 
					else:
						t += '[' + str(k) + ']' + '\t      ' 
				else:
					g = MAX - len(str(int(CD[i][j])))
					t = t + '  ' + str(0)*g + '%.4f'%(CD[i][j])
				n = n +1
			t = t + '\n'
			k += 1
		f.write(t)
		f.close()
	if ((ext == 'P')|(ext == 'A')):
		for i in CD:
			print i
		for i in range(len(D)):
			for j in range(i):
				D[j][i] += D[i][j]
		for i in range(len(D)):
			for j in range(i):
				CD[j][i] += CD[i][j]
		distmat = 'uncorrected_distmat.phy'
		corr_dm = 'corrected_distmat.phy'
		MAX = 0
		for i in CD:
			if (len(str(int(max(i)))) > MAX):
				MAX = len(str(int(max(i))))
		f = open(distmat,'w')
		t = '\t   ' + str(s) + '\n'
		k = 0
		for i in D:
			n = 0
			for j in i:
				if n == 0:
					t += names[k] +  '\t\t  ' + ' '*(2 -len(str(j))) + '%s'%(j)
				else:
					t = t + '\t  ' + ' '*(2 -len(str(j))) + '%s'%(j)
				n = n +1
			t = t + '\n'
			k +=1
		f.write(t)
		f.close()

		f = open(corr_dm,'w')
		t = '\t   ' + str(s) + '\n'
		k = 0
		for i in CD:
			n = 0
			for j in i:
				if n == 0:
					g = MAX - len(str(int(j)))
					t += names[k] +  '\t\t  ' + str(0)*g + '%.4f'%(j)
				else:
					g = MAX - len(str(int(j)))
					t = t + '   ' + str(0)*g + '%.4f'%(j)
				n = n +1
			t = t + '\n'
			k +=1
		f.write(t)
		f.close()
	if ((ext == 'N')|(ext == 'A')):
		D_ut = [[D[i][j] for i in range(len(D))] for j in range(len(D))]
		CD_ut = [[CD[i][j] for i in range(len(CD))] for j in range(len(CD))]
		distmat = 'uncorrected_distmat.nex'
		corr_dm = 'corrected_distmat.nex'
		f = open(distmat,'w')
		t = ('BEGIN DISTANCES;\n DIMENSIONS\n  ntax=' + str(s) + '\n\n ;\n ' +
			'FORMAT\n  triange=UPPER\n\n ;\n Matrix\n')
		k = 1
		for i in range(len(D_ut)):
			n = 0
			for j in range(i,len(D_ut)):
				if n == 0: 
					t += '  '+ 'taxon_' + str(k) + ' '*(1 + -1*(len(str(k))-1)) + '   '*(i)  
				else:
					t += ' ' +' '*(2-len(str(D_ut[i][j]))) + '%s'% (D_ut[i][j])
				n = n +1
			t = t + '\n'
			k +=1
		t += '\n  ;\nEND;'
		f.write(t)
		f.close()

		f = open(corr_dm,'w')
		t = ('BEGIN DISTANCES;\n DIMENSIONS\n  ntax=' + str(s) + '\n\n ;\n ' +
			'FORMAT\n  triange=UPPER\n\n ;\n Matrix\n')
		MAX = 0
		for i in CD:
			if (len(str(int(max(i)))) > MAX):
				MAX = len(str(int(max(i))))
		k = 1
		for i in range(len(CD_ut)):
			n = 0
			for j in range(i,len(CD_ut)):
				if n == 0:
					g = MAX - len(str(int(CD_ut[i][j])))
					t += '  '+ 'taxon_' + str(k) + ' '*(1 + -1*(len(str(k))-1)) + '         '*(i) 
				else:
					g = MAX - len(str(int(CD_ut[i][j])))
					t +=  '  ' + str(0)*g + '%.4f'% (CD_ut[i][j])
				n = n +1
			t = t + '\n'
			k +=1
		t += '\n  ;\nEND;'
		f.write(t)
		f.close()

	

	for i in range(len(admats)):
		fw = open(admats[i],'w')
		t = ""
		k = 0
		for j in M[i]:
			n = 0
			for b in j:
				if n == 0:
					e = 10 - len(rseq[k])
					t = t + rseq[k] + '%*d'%(e,b)
				else:
					t = t + '%*s'%(3,b)
				n = n + 1
			t = t + "\n"
			k += 1
		fw.write(t)
		fw.close()

	
	f = open('fplus_admat.txt','w')
	t = ""
	k = 0
	for i in fp:
		n = 0
		for j in i:
			if n == 0:
				t = t + rseq[k] + '             ' + '%*s'%(8,j)
			else:
				t = t + '   %*s'%(6,j)
			n = n +1
		t = t + '\n'
		k += 1
	f.write(t)
	f.close()

	f = open('fneg_admat.txt','w')
	t = ""
	k = 0
	for i in fn:
		n = 0
		for j in i:
			if n == 0:
				t = t + rseq[k] + '             ' + '%*s'%(8,j)
			else:
				t = t + '   %*s'%(6,j)
			n = n +1
		t = t + '\n'
		k += 1
	f.write(t)
	f.close()
	
	