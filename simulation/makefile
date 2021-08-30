compiler=gfortran
options=$(openmp_threads) -O3 -fmax-stack-var-size=0 -fbounds-check #-heap-arrays
check=-mkl=sequential -check all -warn all -gen_interfaces -fpe0 -ftrapuv -traceback -g
ifort_check=-mkl=parallel -check all -warn all -gen_interfaces -fpe0 -ftrapuv -traceback -g
gfortran_check=-Wall -fbounds-check -ffpe-trap=invalid,overflow -fbacktrace -g
fftw=-lfftw3 -lm -I ./fftw/include -L ./fftw/lib # Single-threaded FFTW
fftw_threads=-lfftw3_threads -lfftw3 -lm -lpthread # Standalone without OpenMP
openmp_threads=-lfftw3_omp -lfftw3 -lm -lpthread -fopenmp -I ./fftw/include -L ./fftw/lib # With OpenMP
all:mathf.mod io.mod setting.mod gpe.mod a.out constants.mod

constants.mod constants.o: constants.f90
	$(compiler) -c constants.f90 $(options)

mathf.mod mathf.o: mathf.f90 constants.mod
	$(compiler) -c mathf.f90 $(options)

io.mod io.o: io.f90 constants.mod
	$(compiler) -c io.f90 $(options)

setting.mod setting.o: setting.f90 constants.mod mathf.mod
	$(compiler) -c setting.f90 $(options)

gpe.mod gpe.o: gpe.f90 constants.mod
	$(compiler) -c gpe.f90 $(options)

a.out: gpe.o io.o setting.o mathf.o
	$(compiler) -o a.out *.o $(options)
	rm *.o *.mod
clean:
	rm *.o *.mod a.out

debug.o: debug.f90 io.o setting.o mathf.o constants.o
	$(compiler) -c debug.f90 $(options)
debug: debug.o
	$(compiler) -o debug.out *.o $(options)
