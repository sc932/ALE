DEFAULT: all

.PHONY: all clean test

all:
	 make -C src

clean:
	make -C src clean

test:
	make -C src test
