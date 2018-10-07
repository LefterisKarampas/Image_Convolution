all:
	mpicc -o build/conv src/main.c src/halo_points.c -lm


mpi:
	mpicc -o build/conv src/main.c src/halo_points.c -lm
	./askdate.sh 4  > /dev/null

openmp:
	mpicc -o build/conv src/main.c src/halo_points.c -lm -fopenmp
	./askdate.sh 1 > /dev/null

clean:
	rm -rf ./build/*

.PHONY: all mpi openmp clean