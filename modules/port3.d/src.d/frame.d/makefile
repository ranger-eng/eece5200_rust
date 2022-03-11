F = gfortran
C = -c -std=legacy
LIB= libport3.a
OBJS = A9RNTC.o  D1MACH.o  EPRINT.o  FRMATR.o  ICEIL.o   ISTKMD.o  NERROR.o  S88FMT.o \
A9RNTD.o  DCEIL.o   ERROFF.o  I0TK00.o  IDCEIL.o  ISTKQU.o  R1MACH.o  SDUMP.o \
A9RNTI.o  DFLR.o    FDUMP.o   I10WID.o  IDFLR.o   ISTKRL.o  RETSRC.o  SETERR.o \
A9RNTL.o  E9RINT.o  FLR.o     I1MACH.o  IFLR.o    ISTKST.o  S1MACH.o  STKDMP.o \
A9RNTR.o  ENTER.o   FRMATD.o  I8SAVE.o  ISTKGT.o  LEAVE.o   S2MACH.o  U9DMP.o \
CEIL.o    ENTSRC.o  FRMATI.o  I8TSEL.o  ISTKIN.o  N5ERR.o   S3MACH.o

.SUFFIXES: .f 

frame: $(LIB) 
	touch frame

$(LIB): $(OBJS)
#        lib /out:../$(LIB)  $(OBJS)  
	ar r ../libport3.a $(OBJS)

clobber:; rm *.o

.f.a:
	$F $C $<
	ar rv ../libport3.a  ?*.o 

.f.o:
	$F $C $<

