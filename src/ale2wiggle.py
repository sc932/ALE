#!/usr/bin/python

def convertToWiggle(inFile):
	kwig = open(inFile + "-kmer.wig", 'w')
	pwig = open(inFile + "-place.wig", 'w')
	iwig = open(inFile + "-insert.wig", 'w')
	dwig = open(inFile + "-depth.wig", 'w')
	wig = open(inFile + ".wig", 'w')

	for line in file(inFile):
		line = line.rstrip()
		if line[0] == "#":
			sp = line.split()
			if sp[1] == "Reference:":
				pwig.write("track name=ALE-place color=0,0,255 group=ALE priority=1\nfixedStep chrom=" + sp[2] + " start=0 step=1\n")
				iwig.write("track name=ALE-insert color=0,0,255 group=ALE priority=1\nfixedStep chrom=" + sp[2] + " start=0 step=1\n")
				dwig.write("track name=ALE-depth color=255,0,0 group=ALE priority=2\nfixedStep chrom=" + sp[2] + " start=0 step=1\n")
				kwig.write("track name=ALE-kmer color=0,255,0 group=ALE priority=3\nfixedStep chrom=" + sp[2] + " start=0 step=1\n")
				wig.write("track name=depth color=0,0,0\nfixedStep chrom=" + sp[2] + " start=0 step=1\n")
			continue
		contig,position,depth,depthLike,placeLike,insertLike,kmerLike = line.split()
		kwig.write(str(kmerLike) + "\n")	
		pwig.write(str(placeLike) + "\n")
		iwig.write(str(insertLike) + "\n")	
		dwig.write(str(depthLike) + "\n")	
		wig.write(str(depth) + "\n")

def main(argv):
        """ """
        convertToWiggle(argv[1])

if __name__ == "__main__":
    import sys
    import getopt
    sys.exit(main(sys.argv))
    

