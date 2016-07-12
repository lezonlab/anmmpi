CC=gcc
PCC=mpicc
FC=gfortran
STRICT=-Wall -Wstrict-prototypes -ansi -pedantic 

anmmpi : anmmpi.o nrstuff.o blzsolvesparse.o MA47.o
	$(PCC) $^ -o $@ -lm -llapack -lblas -lblzpack -lgfortran 

anmmpi.o : anmmpi.c
	$(PCC) -c $< -o $@ -O3

blzsolvesparse.o : blzsolvesparse.c  
	$(CC) -c $< -o $@ -funroll-all-loops -O3

MA47.o : MA47.f
	$(FC) -c $< -o $@ -O3

%.o : %.c
	$(CC) -c $< -o $@ -O3

.PHONY : clean
clean :
	-rm *.o

