from mgenebloomfilter import BloomFilter
import numpy as np
import cProfile as c

bases = {0:"T",1:"A",2:"G",3:"C"}
def genrandomseq(length,total_n):
	cnt = 0
	while cnt < total_n:
		cnt +=1
		yield "".join([bases[i] for i in np.random.random_integers(low=0,high=3,size=length)])



# Testing params
k = 32
n_to_add = 100000
kmers = [i for i in genrandomseq(k,n_to_add)]

n_blooms = 10

def old_f():
	blooms = {}
	for b in xrange(n_blooms):
		blooms[b] = BloomFilter(10000000, 0.01000)
	
	for kmer in kmers:
		_ = [blooms[b].add(kmer) for b in blooms]



def split_f():
	blooms = {}
	for b in xrange(n_blooms):
		blooms[b] = BloomFilter(10000000, 0.01000)
	
	for kmer in kmers:
		hashes = blooms[0].get_all_hashes(kmer)
		_ = [blooms[b].add_from_all_hashes(hashes) for b in blooms]


c.run("old_f(k,n_to_add)")
c.run("split_f(k,n_to_add)")


