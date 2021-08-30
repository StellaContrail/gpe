compiler=gfortran
options=$(fftw) $(silo) -O3 -fmax-stack-var-size=0 -fbounds-check
fftw=-lfftw3 -lm -I ../gpe_fftw/fftw/include -L ../gpe_fftw/fftw/lib
silo=-lsilo -lm -I/usr/local/include
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
