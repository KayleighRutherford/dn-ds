from fractions import Fraction

syn_list = [] #list of synonymous sites per codon
nonsyn_list=[] #list of non synonymous sites per codon
BASES = {'A', 'G', 'T', 'C'}
ns=0 #total non synonymous sites
ss=0 #synonymous sites

codontoaa = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
	'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
	} #dictionary of codons to amino acid translations, with codons as key

def synsites(codon):
	for i in range(len(codon)):
		base = codon[i]
		other_bases = BASES - {base} #all other bases, other than the current one
		syn = 0
		nonsyn=0
		for new_base in other_bases:
			new_codon = codon[:i] + new_base + codon[i + 1:] #replace current base with every other possibility
			if codontoaa[new_codon]==codontoaa[codon]: #if it pulls the same amino acid as before
				syn+=1 #count as synonymous site
			else:
				nonsyn+=1 #otherwise it is non synonymous
		syn_list.append(Fraction(syn, 3)) 
		nonsyn_list.append(Fraction(nonsyn,3))
	return syn_list,nonsyn_list

def dnds(seq1ref,seq2):
	s1=seq1ref #the reference sequence, against which the other sequence is compared
	s2=seq2
	dn=0
	ds=0
	
	for i in xrange(0,len(seq1ref)-3,3): #read sequence in steps of 3 (codons)
		frame1=s1[i:i+3] 
		frame2=s2[i:i+3]
		
		ss,ns=synsites(frame1) #count the number of synonymous and non synonymous sites in the reference sequence codon
		
		
		if frame1 != frame2: #codons are different so there has been a mutation
			for j in xrange(0,3):
				base1=frame1[j] #check each base
				base2=frame2[j] 
				if base1!=base2: #if bases are different
					testmutant=list(frame1)
					testmutant[j]=base2 #replace the reference sequence base with the mutant base
					if codontoaa["".join(testmutant)]==codontoaa[frame1]: #if it pulls the same amino acid as it originally did, the mutation is synonymous
						ds+=1
					else:
						dn+=1
				else:
					continue
		else: #codons are the same so there are no mutations
			continue
	ss=sum(ss)
	ns=sum(ns)		
	print float((dn/ns)/(ds/ss))

			
			
dnds("ATGTCGTTTACTTTGACCAACAAGAACGTGATTTTCGTTGCCGGTCTGGGAGGCATTGGTCTGGACACCAGCAAGGAGCTGCTCAAGCGCGATCTGAAGAACCTGGTGATCCTCGACCGCATTGAGAACCCGGCTGCCATTGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACCTTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGTTGCTGAAGACCATCTTCGCCCAGCTGAAGACCGTCGATGTCCTGATCAACGGAGCTGGTATCCTGGACGATCACCAGATCGAGCGCACCATTGCCGTCAACTACACTGGCCTGGTCAACACCACGACGGCCATTCTGGACTTCTGGGACAAGCGCAAGGGTGGTCCCGGTGGTATCATCTGCAACATTGGATCCGTCACTGGATTCAATGCCATATACCAGGTGCCCGTCTACTCCGGCACCAAGGCCGCCGTGGTCAACTTCACCAGCTCCCTGGCGAAACTGGCACCCATCACCGGCGTGACCGCTTACACCGTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGTTGGATGTTGAGCCCCAGGTTGCTGAGAAGCTCCTGGCTCATCCCACCCAGCCATCGTTGGCCTGCGCCGAGAACTTCGTCAAGGCTATCGAGCTGAACCAGAACGGAGCCATCTGGAAACTGGACTTGGGCAACCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATCTAA","ATGGCGTTTACCTTGACCAACAAGAACGTGGTTTTCGTGGCCGGTCTGGGAGGCATTGGTCTGGACACCAGCAAGGAGCTGGTCAAGCGGGACCTGAAGAACCTGGTGATCCTCGACCGCATTGAGAACCCGGCTGCCATCGCCGAGCTGAAGGCAATCAATCCAAAGGTGACCGTCACCTTCTACCCCTATGATGTGACCGTGCCCATTGCCGAGACCACCAAGCTGCTGAAGACCATCTTCGCCCAGCTGAAGACCATCGATGTCCTGATCAACGGAGCTGGCATCCTGGACGATCACCAGATCGAGCGCACCATCGCCGTCAACTACACCGGCCTGGTGAACACCACGACTGCCATCCTGGACTTCTGGGACAAGCGCAAGGGTGGACCCGGTGGTATCATCTGCAACATTGGATCCGTGACTGGATTCAACGCCATCTACCAGGTGCCCGTTTACTCCGGCACCAAGGCTGCCGTGGTCAACTTCACCAGCTCCCTGGCGAAACTGGCCCCCATCACCGGCGTGACCGCTTACACCGTGAACCCCGGCATCACCCGCACCACCCTGGTGCACAAGTTCAACTCCTGGCTGGATGTTGAGCCCCAGGTGGCCGAGAAGCTCCTGGCTCACCCCACCCAGCCCTCGTTGGCCTGCGCCCAGAACTTTGTGAAGGCCATCGAGCTGAACCAGAACGGTGCCATCTGGAAACTGGACTTGGGCATCCTGGAGGCCATCCAGTGGACCAAGCACTGGGACTCCGGCATCTAA")
		
with open("mel-all.fasta.txt") as m: #read sequence files into lists
	mseq=m.read().splitlines()
m.close()

with open("sim-all.fasta.txt") as s:
	sseq=s.read().splitlines()
s.close()

with open("yak-all.fasta.txt") as y:
	yseq=y.read().splitlines()
y.close()

total=0
for a in mseq: #for every individial in a species list, compare to every individual in another species list. Do for all combinations of species
	for b in sseq:
		total+=dnds(a,b)
print total/72

total=0
for j in sseq:
	for k in yseq:
		total+=dnds(j,k)
print total/72
		
total=0
for x in mseq:
	for z in yseq:
		total+=dnds(x,z)	
print total/144
		
		
		
	