CC = mpicc
OPTIONS = -O3
CFLAGS = -fopenmp
LDFLAGS = -fopenmp
LDLIBS =

main.exe: main.o matmatthread.o c_timer.o
	$(CC) $(OPTIONS) -o main.exe main.o matmatthread.o c_timer.o $(LDFLAGS) $(LDLIBS)

main.o: main.c c_timer.c
	$(CC) $(OPTIONS) -c main.c $(CFLAGS)

matmatthread.o: matmatthread.c
	$(CC) $(OPTIONS) -c matmatthread.c $(CFLAGS)

c_timer.o: c_timer.c
	$(CC) $(OPTIONS) -c c_timer.c $(CFLAGS)

clean:
	rm -f *.o *.exe