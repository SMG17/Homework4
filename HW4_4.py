import pysam
from scipy.special import binom

def fasta_read(genome_file):
	names_line_num = []
	names = []
	genome = []
	for i, line in enumerate(genome_file):
		line = line.strip()
		genome.append(line)
		if line.startswith(">"):
			names_line_num.append(i)
			chromosome = line[1:]
			names.append(chromosome)
	names_line_num.append(len(genome))
	sequences = []
	for i in range(len(names_line_num)-1):
		seq = [''.join(genome[names_line_num[i]+1:names_line_num[i+1]])]
		sequences.extend(seq)
	dictionary = dict(zip(names, sequences))
	return dictionary

def genotype_likelihood(reference_nb, alternate_nb, e=0.01):
	n = alternate_nb + reference_nb
	gen_likelihood_reference = binom(n, alternate_nb)*e**alternate_nb*(1.0-e)**reference_nb
	gen_likelihood_alternate = binom(n, alternate_nb)*(1.0-e)**alternate_nb*e**reference_nb
	return gen_likelihood_reference, gen_likelihhod_alternate

yeast_chr = open("yeast_chr.fa")
yeast_genome = fasta_read(yeast_chr)

bamfile = pysam.AlignmentFile("output.sorted.bam")

output = open("SNP_genotyoe_likelihoods","w")

for pileupcolumn in bamfile.pileup():
	chromosome = pileupcolumn.reference_name
	position = pileupcolumn.pos
	reference_allele = yeast_genome[chromosome][position]
	if pileupcolumn.n > 100: continue
	reference_nb = 0
	alternate_alleles = {}
	for pileupread in pileupcolumn.pileups:
		if pileupread.is_del or pileupread.is_refskip: continue
		read = pileupread.alignment.query_sequence[pileupread.query_position]
		if read != reference_allele:
			if read in alternate_alleles:
				alternate_alleles[read] += 1
			else:
				alternate_alleles[read] = 1
		else:
			reference_nb += 1
	alternate_nb = 0
	alternate_allele = ""
	for allele in alternate_alleles:
		if alternate_alleles[allele] > alternate_nb:
			alternate_nb = alternate_alleles[allele]
			alternate_allele = allele
	if alternate_nb > 0:
		gen_likelihood_ref, gen_likelihood_alt = genotype_likelihood(reference_nb, alternate_nb)
		output.write("%s\t%d\t%s\t%s\t%f\t%f\n" % (chromosome, position, reference_allele, alternate_allele, gen_likelihood_ref, gen_likelihood_alt))		

yeast_chr.close()
bamfile.close()
output.close()
				

