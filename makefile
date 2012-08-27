DEFAULT: all

all:
	 make -C src

clean:
	make -C src clean

test:
	make -C src test
