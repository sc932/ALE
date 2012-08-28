SAMTOOLS_PATH := samtools-0.1.18
 
UNAME := $(shell uname)

CC := gcc

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

.PHONY: all clean cleantest test

samlib:
	make -C $(SAMTOOLS_PATH) lib

$(SAMTOOLS_PATH)/bam.o : samlib
		true

ALElike.o: ALElike.c ALElike.h
		$(CC) -c $(CFLAGS) $(GNU_FLAG) ALElike.c -I$(SAMTOOLS_PATH)

ALEhelpers.o: ALEhelpers.c ALEhelpers.h 
		$(CC) -c $(CFLAGS) $(GNU_FLAG) ALEhelpers.c -I$(SAMTOOLS_PATH)

ALE: ALE.c ALE.h ALEhelpers.o ALElike.o $(SAMTOOLS_PATH)/bam.o
		$(CC) $(CFLAGS) $(GNU_FLAG) ALE.c ALElike.o ALEhelpers.o -o ALE -lz -lm -I$(SAMTOOLS_PATH) -L$(SAMTOOLS_PATH) -lbam

synthReadGen: synthReadGen.c
		$(CC) $(CFLAGS) synthReadGen.c -o synthReadGen -lz -lm -I$(SAMTOOLS_PATH) -L$(SAMTOOLS_PATH)

readFileSplitter: readFileSplitter.c
		$(CC) $(CFLAGS) readFileSplitter.c -o readFileSplitter

GCcompFinder: GCcompFinder.c
		$(CC) $(CFLAGS) GCcompFinder.c -o GCcompFinder

all: GCcompFinder readFileSplitter synthReadGen ALE

install: ALE synthReadGen

clean: cleantest
		rm -f *.o GCcompFinder readFileSplitter synthReadGen2 ALE
		make -C $(SAMTOOLS_PATH) clean

cleantest: 
		rm -f ../test/*-testPhiX.ale* ../test/*-tinytestPhiX.ale* ../test/*-tinytestPhiX-broken.ale*

test: ALE cleantest test1 test2 test3 test4 test-single-relations test-pair-relations test-pair-chimer-baseline test-pair-chimer

%-testPhiX.ale: % ALE
		./ALE $< ../test/PhiX174.fsa $@ > $@.log 2>&1 && ls -l $@

%-tinytestPhiX.ale: % ALE ../test/PhiX-tiny.sam-testPhiX.ale
		./ALE -pm ../test/PhiX-tiny.sam-testPhiX.ale.param $< ../test/PhiX174.fsa $@ > $@.log 2>&1 && ls -l $@

%-tinytestPhiX-broken.ale: % ALE ../test/PhiX-tiny.sam-testPhiX.ale
		./ALE -pm ../test/PhiX-tiny.sam-testPhiX.ale.param $< ../test/PhiX174-broken.fsa $@ > $@.log 2>&1 && ls -l $@

test1: ../test/PhiX-sort.bam-testPhiX.ale
		diff -q $< ../test/PhiX.ale > $<.diff && echo "$@ PASSED"

test2: ../test/PhiX.bam-testPhiX.ale
		diff -q $< ../test/PhiX.ale > $<.diff && echo "$@ PASSED"

test3: ../test/PhiX.sam-testPhiX.ale
		diff -q $< ../test/PhiX.ale > $<.diff && echo "$@ PASSED"

test4: ../test/PhiX-tiny.sam-testPhiX.ale
		diff -q $< ../test/PhiX-tiny.ale > $<.diff && echo "$@ PASSED"

test-single-relations: ../test/PhiX-tiny.sam-testPhiX.ale  ../test/PhiX-tiny-partial-single.sam-testPhiX.ale  ../test/PhiX-tiny-partial2-single.sam-testPhiX.ale  ../test/PhiX-tiny-partial3-single.sam-testPhiX.ale  ../test/PhiX-tiny-partial4-single.sam-testPhiX.ale  ../test/PhiX-tiny-partial5-single.sam-testPhiX.ale  ../test/PhiX-tiny-+unmapped-single.sam-testPhiX.ale  ../test/PhiX-tiny-unmapped-single.sam-testPhiX.ale
		awk '/ALE_score/ {print}' $^ > /tmp/test-relations.alescore && sort -k3 -n -r /tmp/test-relations.alescore | diff - /tmp/test-relations.alescore && rm -f /tmp/test-relations.alescore && echo "$@ PASSED"
		
test-pair-relations: ../test/PhiX-tiny-stretch0-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny.sam-testPhiX.ale ../test/PhiX-tiny-stretch1-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch2-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch3-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch4-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch5-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch6-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch7-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch8-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-removed-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-removed2-pair.sam-tinytestPhiX.ale
		awk '/ALE_score/ {print}' $^ > /tmp/test-relations.alescore && sort -k3 -n -r /tmp/test-relations.alescore | diff - /tmp/test-relations.alescore && rm -f /tmp/test-relations.alescore && echo "$@ PASSED"

test-pair-chimer-baseline: ../test/PhiXbroken-tiny-stretch0-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch1-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch2-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch3-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch4-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch5-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch6-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch7-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch8-pair.sam-tinytestPhiX-broken.ale
		awk '/ALE_score/ {print}' $^ > /tmp/test-relations.alescore && sort -k3 -n -r /tmp/test-relations.alescore | diff - /tmp/test-relations.alescore && rm -f /tmp/test-relations.alescore && echo "$@ PASSED"

test-pair-chimer: ../test/PhiXbroken-tiny-stretch0-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch1-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch2-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch3-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch4-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch5-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch6-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-chimer-pair.sam-tinytestPhiX-broken.ale # ../test/PhiXbroken-tiny-removed-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-removed2-pair.sam-tinytestPhiX-broken.ale 
		awk '/ALE_score/ {print}' $^ > /tmp/test-relations.alescore && sort -k3 -n -r /tmp/test-relations.alescore | diff - /tmp/test-relations.alescore && rm -f /tmp/test-relations.alescore && echo "$@ PASSED"

