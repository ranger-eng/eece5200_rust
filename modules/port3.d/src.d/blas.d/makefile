F = gfortran
C = -c -std=legacy
LIB= libport3.a
OBJS = CABS1.o SROT2.o CROTG.o   DASUM.o   DROTG.o   IDAMAX.o  ISAMAX.o  SAXPY.o   \
CAMAX.o   CSCAL.o   DAXPY.o   DSCAL.o   IDAMIN.o  ISAMIN.o  SCASUM.o  SROTG.o \
CAXPY.o   CSROT.o   DCOPY.o   DSUM.o    IDMAX.o   ISMAX.o   SCOPY.o   SSCAL.o \
CCOPY.o   CSSCAL.o  DDOT.o    DSWAP.o   IDMIN.o   ISMIN.o   SDOT.o    SSUM.o \
CDOTC.o   CSWAP.o   DNRM2.o   ICAMAX.o  IIMAX.o   SAMAX.o   SNRM2.o   SSWAP.o \
CDOTU.o   DAMAX.o   DROT.o    ICAMIN.o  IIMIN.o   SASUM.o   SROT.o

.SUFFIXES: .f 

blas: $(LIB) 
	touch blas

$(LIB): $(OBJS)
#       lib /out:../$(LIB) /verbose $(OBJS)
	ar r ../libport3.a $(OBJS)

clobber:; rm *.o

.f.a:
	$F $C $<
	ar r ../libport3.a $*.o 

.f.o:
	$F $C $<

