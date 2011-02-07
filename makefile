install: ALE synthReadGen2

ALE: ALE.c ALE.h geneTree.h
		cc -g -O2 ALE.c -o ALE -lz -lm

synthReadGen2: synthReadGen2.c
		cc -g -O2 synthReadGen2.c -o synthReadGen2 -lz -lm