LEV1 = frame blas GESS EIGEN util
MAKE =  /usr/bin/make
LIB =	libport3.a

all:;   $(MAKE)  $(LEV1)

$(LEV1):;
	-cd src.d/$@.d; $(MAKE) $@
lib:

install:;

clean:
	rm -f  src.d/blas.d/*.o
	rm -f  src.d/frame.d/*.o
	rm -f  src.d/EIGEN.d/*.o
	rm -f  src.d/GESS.d/*.o
	rm -f  src.d/util.d/*.o
