CC = mpicc
OPTIONS = -O3
CFLAGS = -fopenmp
LDFLAGS = -fopenmp
LDLIBS =

main.out: mainmatmatdist.o matmatdist.o c_timer.o
	$(CC) $(OPTIONS) -o main.out mainmatmatdist.o matmatdist.o c_timer.o $(LDFLAGS) $(LDLIBS)

mainmatmatdist.o: mainmatmatdist.c
	$(CC) $(OPTIONS) -c mainmatmatdist.c $(CFLAGS)

matmatdist.o: matmatdist.c
	$(CC) $(OPTIONS) -c matmatdist.c $(CFLAGS)

c_timer.o: c_timer.c
	$(CC) $(OPTIONS) -c c_timer.c $(CFLAGS)

clean:
	rm -f *.exe *.out *.o