def fasta_read(x):
	names_line_num = []
	names = []
	genome = []
	for i, line in enumerate(x):
		line = line.strip()
		genome.append(line)
		if line.startswith(">"):
			names_line_num.append(i)
			chromosome = line[1:]
			chromosome = chromosome.replace("_"," ")
			names.append(chromosome)

	names_line_num.append(len(genome))
	sequences = []
	for i in range(len(names_line_num)-1):
		seq = [''.join(genome[names_line_num[i]+1:names_line_num[i+1]])]
		sequences.extend(seq)

	dictionary = dict(zip(names, sequences))
	return dictionary

x = open("/mnt/c/Users/SMG/Desktop/my_genome_ex.txt")
genome_dictionary = fasta_read(x)

# 1-a)
sequences = []
BWT_dic = {}
for key, value in genome_dictionary.items():
	seq = value + "$"
	sequences.append(seq)
	new_sequences = []
	for i in range(len(seq)):
		new_seq = seq[len(seq)-1] + seq[0:(len(seq)-1)]		
		seq = new_seq
		new_sequences.append(new_seq)
		BWT_dic[key] = new_sequences

BWT_sorted_dic = {}
for key, value in BWT_dic.items():
	BWT_sorted_dic[key] = sorted(value)

BWT_last_column = {}
for key, value in BWT_sorted_dic.items():
	last_characters = []
	for i in range(len(value)):
		last_characters.extend(BWT_sorted_dic[key][i][-1])
	BWT_last_column[key] = ''.join(last_characters)

for key, value in BWT_last_column.items():
	myOutput = open("BWT_%s.txt" % key, "w")
	myOutput.write("The BWT of %s is %s." % (key, value))
	myOutput.close()

# 1-b)
SA_dic = {}
for key, value in BWT_sorted_dic.items():
	SA = []
	myOutputSA = open("CheckpointedSA_%s.txt" % key,"w")
	for i in range(int(round(float(len(value))/2))):
		suffix = BWT_sorted_dic[key][i*2]
		array = suffix[suffix.find("$")+1:]
		arrayLength = len(array)
		SA.append(arrayLength)
		myOutputSA.write("%s	%s\n" % (i*2, arrayLength))
	myOutputSA.close()
	SA_dic[key] = SA

# 1-c)
for key, value in BWT_last_column.items():
	value_sorted = sorted(value)
	value_sorted = ''.join(value_sorted)
	dollar_Cc = len(value_sorted[:value_sorted.find("$")])
	A_Cc = len(value_sorted[:value_sorted.find("A")])
	C_Cc = len(value_sorted[:value_sorted.find("C")])
	G_Cc = len(value_sorted[:value_sorted.find("G")])
	T_Cc = len(value_sorted[:value_sorted.find("T")])
	myOutputCc = open("TableC[c]_%s.txt" % key,"w")
	myOutputCc.write("$	A	C	G	T\n%s	%s	%s	%s	%s" % (dollar_Cc, A_Cc, C_Cc, G_Cc, T_Cc))
	myOutputCc.close()

# 1-d)
for key, value in BWT_last_column.items():
	dollar_Occ = []
	dollar = 0
	A_Occ = []
	A = 0
	C_Occ = []
	C = 0
	G_Occ = []
	G = 0
	T_Occ = []
	T = 0
	for i in range(len(value)):
		if value[i] == '$':
			dollar += 1
			dollar_Occ.append(dollar)
		else:
			dollar_Occ.append(dollar)
	for i in range(len(value)):
		if value[i] == 'A':
			A += 1
			A_Occ.append(A)
		else:
			A_Occ.append(A)
	for i in range(len(value)):
		if value[i] == 'C':
			C += 1
			C_Occ.append(C)
		else:
			C_Occ.append(C)
	for i in range(len(value)):
		if value[i] == 'G':
			G += 1
			G_Occ.append(G)
		else:
			G_Occ.append(G)
	for i in range(len(value)):
		if value[i] == 'T':
			T += 1
			T_Occ.append(T)
		else:
			T_Occ.append(T)
	myOutputOcc = open("OccTable_%s.txt" % key, "w")
	for i in range(int(round(float(len(value))/2))):
		myOutputOcc.write("%s	%s	%s	%s	%s" % (i*2, A_Occ[i*2], C_Occ[i*2], G_Occ[i*2], T_Occ[i*2]))
	myOutputOcc.close()

x.close()
