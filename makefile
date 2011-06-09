SamtoolsPath = /jgi/tools/misc_bio/samtools/DEFAULT

install: ALE synthReadGen2

ALE: ALEv2.c ALE.h geneTree.h ALEv2.h ALElike.h
		cc -g -O2 ALEv2.c -o ALE -lz -lm -I$(SamtoolsPath)/include

synthReadGen2: synthReadGen2.c
		cc -g -O2 synthReadGen2.c -o synthReadGen2 -lz -lm

readFileSplitter: readFileSplitter.c
		cc -g -O2 readFileSplitter.c -o readFileSplitter

GCcompFinder: GCcompFinder.c
		cc -g -O2 GCcompFinder.c -o GCcompFinder

all: GCcompFinder readFileSplitter synthReadGen2 ALE

DEFAULT: all

clean:
		rm -f GCcompFinder readFileSplitter synthReadGen2 ALE


		