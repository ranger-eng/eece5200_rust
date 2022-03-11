F = gfortran
C = -c -std=legacy
LIB= libport3.a
OBJS = CDDIV.o



.SUFFIXES: .f 

util: $(LIB) 
	touch util

$(LIB): $(OBJS)
	ar r ../libport3.a $(OBJS)

clobber:; rm *.o

.f.a:
	$F $C $<
	ar rv ../libport3.a  ?*.o 

.f.o:
	$F $C $<

