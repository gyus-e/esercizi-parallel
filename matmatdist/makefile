CC = gcc
CFLAGS = -O3 -fopenmp
LDFLAGS =

main.exe: main.o matmatthread.o c_timer.o
	$(CC) $(CFLAGS) -o main.exe main.o matmatthread.o c_timer.o $(LDFLAGS)

main.o: main.c c_timer.c matmatthread.h
	$(CC) $(CFLAGS) -c main.c

matmatthread.o: matmatthread.c matmatthread.h
	$(CC) $(CFLAGS) -c matmatthread.c

c_timer.o: c_timer.c
	$(CC) $(CFLAGS) -c c_timer.c

clean:
	rm -f *.o *.exe