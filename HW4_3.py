import pysam
import matplotlib.pyplot as mpl

bamfile = pysam.AlignmentFile("output.sorted.bam")

list_mapping_quality = []
posStrand = 0
negStrand = 0
for read in bamfile.fetch():
	list_mapping_quality.append(read.mapping_quality)
	if read.is_reverse == True:
		negStrand += 1
	else:
		posStrand += 1
lenReads = negStrand + posStrand
negProp = float(negStrand)/lenReads
posProp = float(posStrand)/lenReads

mpl.hist(list_mapping_quality, bins=100)
mpl.savefig("histogram_3c.jpeg")

print "The proportion of reads from the forward strand is %s" % posProp	
print "The proportion of reads from the reverse strand is %s" % negProp
