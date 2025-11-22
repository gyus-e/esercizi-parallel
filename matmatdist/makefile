CC = gcc
CFLAGS = -O3
LDFLAGS =

main.exe: main.o matmatblock.o matmat.o c_timer.o
	$(CC) $(CFLAGS) -o main.exe main.o matmatblock.o matmat.o c_timer.o $(LDFLAGS)

main.o: main.c c_timer.c matmatblock.h matmat.h
	$(CC) $(CFLAGS) -c main.c

matmatblock.o: matmatblock.c matmatblock.h matmat.h
	$(CC) $(CFLAGS) -c matmatblock.c

matmat.o: matmat.c matmat.h
	$(CC) $(CFLAGS) -c matmat.c

c_timer.o: c_timer.c
	$(CC) $(CFLAGS) -c c_timer.c

clean:
	rm -f *.o *.exe