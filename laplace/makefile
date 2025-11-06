CC = mpicc
OPTIONS =
CFLAGS =
LDFLAGS =
LDLIBS =

laplace.out: mainlaplace.o laplace.o c_timer.o
	$(CC) $(OPTIONS) -o laplace.out mainlaplace.o laplace.o c_timer.o $(LDFLAGS) $(LDLIBS)

mainlaplace.o: mainlaplace.c
	$(CC) $(OPTIONS) -c mainlaplace.c $(CFLAGS)

laplace.o: laplace.c
	$(CC) $(OPTIONS) -c laplace.c $(CFLAGS)

c_timer.o: c_timer.c
	$(CC) $(OPTIONS) -c c_timer.c $(CFLAGS)

clean:
	rm -f *.out *.o