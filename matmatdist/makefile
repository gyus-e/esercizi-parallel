CC = gcc
CFLAGS = -O3
LDFLAGS =

all: main

main: main.o matmat.o c_timer.o
	$(CC) $(CFLAGS) -o main main.o matmat.o c_timer.o $(LDFLAGS)

main.o: main.c c_timer.c matmat.h
	$(CC) $(CFLAGS) -c main.c

matmat.o: matmat.c matmat.h
	$(CC) $(CFLAGS) -c matmat.c

c_timer.o: c_timer.c
	$(CC) $(CFLAGS) -c c_timer.c

clean:
	rm -f *.o main