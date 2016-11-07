.PHONY : test


test :
	mpifort test.f90 gfs_state.o -o test -I../build/include -I/usr/lib64/gfortran/modules -L../build/lib -lletkf -lblas -llapack -lnetcdff
	LD_LIBRARY_PATH=../build/lib   mpirun -mca mpi_yield_when_idle 1 -n 5 ./test
#
