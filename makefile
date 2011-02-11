install: ALE synthReadGen2

ALE: ALEv2.c ALE.h geneTree.h ALEv2.h ALElike.h
		cc -g -O2 ALEv2.c -o ALE -lz -lm

synthReadGen2: synthReadGen2.c
		cc -g -O2 synthReadGen2.c -o synthReadGen2 -lz -lm