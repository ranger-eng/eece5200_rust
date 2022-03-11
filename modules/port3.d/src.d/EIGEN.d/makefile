F = gfortran
C = -c -std=legacy
LIB= libport3.a
OBJS = EIGEN.o HQR2.o ORTHE.o ORTRA.o \
	DEIGEN.o DHQR2.o DORTHE.o DORTRA.o



.SUFFIXES: .f 

EIGEN: $(LIB) 
	touch eigen

$(LIB): $(OBJS)
	ar r ../libport3.a $(OBJS)

clobber:; rm *.o

.f.a:
	$F $C $<
	ar rv ../libport3.a  ?*.o 

.f.o:
	$F $C $<

