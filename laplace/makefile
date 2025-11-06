CC = mpicc
CFLAGS = -std=gnu90
LDLIBS = -lm

laplace.out: mainlaplace.c laplace.c c_timer.c
	$(CC) -O3 -o laplace.out mainlaplace.c laplace.c c_timer.c

clean:
	rm -f *.out *.o *.txt *.err