SAMTOOLS_PATH= /jgi/tools/misc_bio/samtools/DEFAULT
CFLAGS= -g -O2


install: ALE synthReadGen2

ALE: ALEv2.c ALE.h geneTree.h ALEv2.h ALElike.h
		$(CC) $(CFLAGS) ALEv2.c -o ALE -lz -lm -I$(SAMTOOLS_PATH)/include -L$(SAMTOOLS_PATH)/lib -lbam -Wall

synthReadGen2: synthReadGen2.c
		$(CC) $(CFLAGS) synthReadGen2.c -o synthReadGen2 -lz -lm

readFileSplitter: readFileSplitter.c
		$(CC) $(CFLAGS) readFileSplitter.c -o readFileSplitter

GCcompFinder: GCcompFinder.c
		$(CC) $(CFLAGS) GCcompFinder.c -o GCcompFinder

all: GCcompFinder readFileSplitter synthReadGen2 ALE

DEFAULT: all

clean:
		rm -f GCcompFinder readFileSplitter synthReadGen2 ALE


		