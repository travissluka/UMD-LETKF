.PHONY : test


test :
	mpif90 test.f90 -o test -I../build/include -L../build/lib -lletkf -lblas -llapack
	LD_LIBRARY_PATH=../build/lib mpirun -n 5 ./test
