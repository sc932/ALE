DEFAULT: all

.PHONY: all clean test example

all:
	 make -C src

clean:
	make -C src clean

test:
	make -C src test

example:
	make -C src example
