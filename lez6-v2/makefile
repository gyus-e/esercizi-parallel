CC = gcc
CFLAGS = -O3 -fopenmp -Wall -Wextra -std=gnu90
LDFLAGS = -O3 -fopenmp
LDLIBS = -lm

a.out: main.o maxsum.o
	$(CC) -o a.out main.o maxsum.o $(LDFLAGS) $(LDLIBS)

main.o: main.c maxsum.h.gch
	$(CC) -c main.c $(CFLAGS)

maxsum.o: maxsum.c maxsum.h.gch
	$(CC) -c maxsum.c $(CFLAGS)

maxsum.h.gch: maxsum.h
	$(CC) -c maxsum.h

clean:
	rm -f *.out *.err *.o *.h.gch
