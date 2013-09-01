CC = gcc
CFLAGS = -O3 -fPIC -march=nocona -ffast-math -fomit-frame-pointer 

#-L/home/joao/matlab/bin/glnx86/
#CC=icc
#CFLAGS = -xP -fast -fPIC
OMPFLAGS = -fopenmp

MEXEXT=a64
MATLABDIR=/usr/local/MATLAB/R2013a
INCLUDES=-I$(MATLABDIR)/extern/include
LDIRS= -L$(MATLABDIR)/bin/glnx$(MEXEXT)

EXE_TARGETS = chi2float chi2double chi2_mex.mex$(MEXEXT)
LIB_TARGETS = libchi2.so
all:	$(EXE_TARGETS) $(LIB_TARGETS)

chi2float:	chi2float.c chi2float.h Makefile
	$(CC) -D__MAIN__  $(CFLAGS) $(OMPFLAGS) -o chi2float chi2float.c

chi2double: chi2double.c chi2double.h Makefile
	$(CC) -D__MAIN__  $(CFLAGS)  $(OMPFLAGS) -o chi2double chi2double.c

libchi2.so:	chi2double.c chi2double.h chi2float.c chi2float.h Makefile
	$(CC) $(CFLAGS) $(OMPFLAGS) -shared -Wl,-soname=libchi2.so -fPIC chi2double.c chi2float.c -o libchi2.so

chi2double.o : chi2double.c chi2double.h Makefile
	$(CC) -D__MAIN__  $(CFLAGS) -c $(OMPFLAGS) -o chi2double.o chi2double.c

chi2float.o: chi2float.c chi2double.h Makefile
	$(CC) -D__MAIN__  $(CFLAGS) -c $(OMPFLAGS) -o chi2float.o chi2float.c

chi2_mex.o: chi2_mex.c
	$(CC) $(CFLAGS) -c $(INCLUDES) $(OMPFLAGS) -o chi2_mex.o chi2_mex.c

chi2_mex.mex$(MEXEXT): 	chi2_mex.c chi2_mex.o chi2double.o
	$(CC) $(OMPFLAGS) chi2_mex.o  $(LDIRS) $(CFLAGS) -lmex -shared -o chi2_mex.mex$(MEXEXT) chi2double.o

#chi2_mex.mex$(MEXEXT):  chi2_mex.c chi2_mex.o chi2float.o
#	$(CC) $(OMPFLAGS) chi2_mex.o  $(LDIRS) $(CFLAGS) -lmex -shared -o chi2_mex.mex$(MEXEXT) chi2float.o

# default installation of libomp cannot be opened using dlopen() as would be required e.g. for Python


clean:
	rm -f *.o $(EXE_TARGETS) $(LIB_TARGETS)

timing:	$(EXE_TARGETS)
	time ./chi2float
	time ./chi2double

