CC = mpicc
CFLAGS = -std=gnu90
LDLIBS = -lm

laplace.out: mainlaplace.c laplace.c c_timer.c laplace_blocking.out
	$(CC) -o laplace.out mainlaplace.c laplace.c c_timer.c

laplace_blocking.out: mainlaplace.c laplace_blocking.c c_timer.c
	$(CC) -o laplace_blocking.out mainlaplace.c laplace_blocking.c c_timer.c

clean:
	rm -f *.out *.o *.txt *.err