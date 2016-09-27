.PHONY : test


test :
	mpifort test.f90 -o test -I../build/include -L../build/lib -lletkf -lblas -llapack -lnetcdff
	LD_LIBRARY_PATH=../build/lib mpirun -n 5 ./test
