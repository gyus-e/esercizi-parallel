CC = mpicc
CFLAGS = -std=gnu90
LDLIBS = -lm

laplace.out: mainlaplace.c laplace.c c_timer.c
	mpicc -o laplace.out mainlaplace.c laplace.c c_timer.c $(CFLAGS) $(LDLIBS)

laplace_blocking.out: mainlaplace.c laplace_blocking.c c_timer.c
	mpicc -o laplace_blocking.out mainlaplace.c laplace_blocking.c c_timer.c $(CFLAGS) $(LDLIBS)

clean:
	rm -f *.out *.o *.txt