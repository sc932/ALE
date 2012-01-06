#!/usr/bin/python

def convertToWiggle(inFile):
	kwig = open(inFile + "-kmer.wig", 'w')
	pwig = open(inFile + "-place.wig", 'w')
	dwig = open(inFile + "-depth.wig", 'w')
	twig = open(inFile + "-total.wig", 'w')
	wig = open(inFile + ".wig", 'w')

	for line in file(inFile):
		line = line.rstrip()
		if line[0] == "#":
			sp = line.split()
			if sp[1] == "Reference:":
				kwig.write("track name=ALE-k\nfixedStep chrom=" + sp[2] + " start=0 step=1\n")
				pwig.write("track name=ALE-p\nfixedStep chrom=" + sp[2] + " start=0 step=1\n")
				dwig.write("track name=ALE-d\nfixedStep chrom=" + sp[2] + " start=0 step=1\n")
				twig.write("track name=ALE-t\nfixedStep chrom=" + sp[2] + " start=0 step=1\n")
				wig.write("track name=depth\nfixedStep chrom=" + sp[2] + " start=0 step=1\n")
			continue
		contig,position,depth,depthLike,placeLike,kmerLike,totalLike = line.split()
		kwig.write(str(kmerLike) + "\n")	
		pwig.write(str(placeLike) + "\n")	
		dwig.write(str(depthLike) + "\n")	
		twig.write(str(totalLike) + "\n")	
		wig.write(str(depth) + "\n")

def main(argv):
        """ """
        convertToWiggle(argv[1])

if __name__ == "__main__":
    import sys
    import getopt
    sys.exit(main(sys.argv))
    

