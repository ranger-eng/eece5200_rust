F = gfortran
C = -c -std=legacy
LIB= libport3.a
OBJS = G4ECE.o GEBS.o GECE.o GEFS.o GELU.o GESS.o \
DG4ECE.o DGEBS.o DGECE.o DGEFS.o DGELU.o DGESS.o

.SUFFIXES: .f 

GESS: $(LIB) 
	touch GESS

$(LIB): $(OBJS)
	ar r ../libport3.a $(OBJS)

clobber:; rm *.o

.f.a:
	$F $C $<
	ar rv ../libport3.a  ?*.o 

.f.o:
	$F $C $<

