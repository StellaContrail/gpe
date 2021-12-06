compiler=gfortran
options=$(fftw) $(silo) -O0 -lm -I/usr/local/include $(check) 
fftw=-lfftw3
silo=-lsilo
check=-g -fcheck=all -Wall -fbounds-check
max_memory=-fmax-stack-var-size=0
all:mathf.mod silo_converter.mod a.out

mathf.mod mathf.o: mathf.f90
	$(compiler) -c mathf.f90 $(options)

silo_converter.mod silo_converter.o: silo_converter.f90
	$(compiler) -c silo_converter.f90 $(options)

a.out: mathf.o silo_converter.o
	$(compiler) -o a.out *.o $(options)
	rm *.o *.mod
clean:
	rm *.o *.mod a.out
