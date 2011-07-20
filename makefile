SAMTOOLS_PATH=samtools-0.1.17

#/jgi/tools/misc_bio/samtools/DEFAULT
CFLAGS= -g -O3

DEFAULT: all

samlib:
	make -C $(SAMTOOLS_PATH) lib

ALE: ALEv2.c ALE.h geneTree.h ALEv2.h ALElike.h samlib
		$(CC) $(CFLAGS) ALEv2.c -o ALE -lz -lm -I$(SAMTOOLS_PATH) -L$(SAMTOOLS_PATH) -lbam -Wall

synthReadGen2: synthReadGen2.c
		$(CC) $(CFLAGS) synthReadGen2.c -o synthReadGen2 -lz -lm

readFileSplitter: readFileSplitter.c
		$(CC) $(CFLAGS) readFileSplitter.c -o readFileSplitter

GCcompFinder: GCcompFinder.c
		$(CC) $(CFLAGS) GCcompFinder.c -o GCcompFinder

all: GCcompFinder readFileSplitter synthReadGen2 ALE

install: ALE synthReadGen2

clean:
		rm -f GCcompFinder readFileSplitter synthReadGen2 ALE


		
