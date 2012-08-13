SAMTOOLS_PATH=samtools-0.1.18

CFLAGS:= -g -O3
# note for Mac OS X and other non-gnu comment out this next variable
GNU_FLAG=-D_GNU_SRC

DEFAULT: all

samlib:
	make -C $(SAMTOOLS_PATH) lib

ALElike.o: ALElike.c ALElike.h
		$(CC) -c $(CFLAGS) $(GNU_FLAG) ALElike.c -I$(SAMTOOLS_PATH)

ALEhelpers.o: ALEhelpers.c ALEhelpers.h 
		$(CC) -c $(CFLAGS) $(GNU_FLAG) ALEhelpers.c -I$(SAMTOOLS_PATH)

ALE: ALE.c ALE.h ALEhelpers.o ALElike.o samlib
		$(CC) $(CFLAGS) $(GNU_FLAG) ALE.c ALElike.o ALEhelpers.o -o ALE -lz -lm -I$(SAMTOOLS_PATH) -L$(SAMTOOLS_PATH) -lbam

synthReadGen: synthReadGen.c
		$(CC) $(CFLAGS) synthReadGen.c -o synthReadGen -lz -lm -I$(SAMTOOLS_PATH) -L$(SAMTOOLS_PATH)

readFileSplitter: readFileSplitter.c
		$(CC) $(CFLAGS) readFileSplitter.c -o readFileSplitter

GCcompFinder: GCcompFinder.c
		$(CC) $(CFLAGS) GCcompFinder.c -o GCcompFinder

all: GCcompFinder readFileSplitter synthReadGen ALE test

install: ALE synthReadGen

clean:
		rm -f *.o GCcompFinder readFileSplitter synthReadGen2 ALE
		make -C $(SAMTOOLS_PATH) clean

test: ALE test1 test2 test3

test1: ALE
		./ALE ../test/PhiX-sort.bam ../test/PhiX174.fsa /tmp/PhiX.ale && diff -q /tmp/PhiX.ale ../test/PhiX.ale && rm -f /tmp/PhiX.ale /tmp/PhiX.ale.param

test2: ALE
		./ALE ../test/PhiX.bam ../test/PhiX174.fsa /tmp/PhiX.ale && diff -q /tmp/PhiX.ale ../test/PhiX.ale && rm -f /tmp/PhiX.ale /tmp/PhiX.ale.param

test3: ALE
		./ALE ../test/PhiX.sam ../test/PhiX174.fsa /tmp/PhiX.ale && diff -q /tmp/PhiX.ale ../test/PhiX.ale && rm -f /tmp/PhiX.ale /tmp/PhiX.ale.param

