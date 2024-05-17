
mpi_rand: mpi_rand.c
	mpicc  -Wall --std=c99  -o $@ $<

clean:
	rm -f mpi_rand	
