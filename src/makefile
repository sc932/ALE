SAMTOOLS_PATH=samtools-0.1.18

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
# do something Linux-y 
GNU_FLAG=-D_GNU_SRC
endif
ifeq ($(UNAME), Solaris)
# do something Solaris-y 
endif
ifeq ($(UNAME), Darwin)
# do something Mac-y
endif

CFLAGS := -g -O3

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

clean: cleantest
		rm -f *.o GCcompFinder readFileSplitter synthReadGen2 ALE
		make -C $(SAMTOOLS_PATH) clean

cleantest: 
		rm -f ../test/*-testPhiX.ale*

test: ALE cleantest test1 test2 test3 test4

%-testPhiX.ale: % ALE
		./ALE $< ../test/PhiX174.fsa $@ && ls -l $@

test1: ../test/PhiX-sort.bam-testPhiX.ale
		diff -q $< ../test/PhiX.ale && echo "$@ PASSED"

test2: ../test/PhiX.bam-testPhiX.ale
		diff -q $< ../test/PhiX.ale && echo "$@ PASSED"

test3: ../test/PhiX.sam-testPhiX.ale
		diff -q $< ../test/PhiX.ale && echo "$@ PASSED"

test4: ../test/PhiX-tiny.sam-testPhiX.ale
		diff -q $< ../test/PhiX-tiny.ale && echo "$@ PASSED"
