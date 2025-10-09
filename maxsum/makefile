CC = gcc
CFLAGS = -O3 -fopenmp -Wall -Wextra -std=gnu90
LDFLAGS = -O3 -fopenmp
LDLIBS = -lm

mainmaxsum.out: mainmaxsum.o maxsum.o c_timer.o
	$(CC) -o mainmaxsum.out mainmaxsum.o maxsum.o c_timer.o $(LDFLAGS) $(LDLIBS)

mainmaxsum.o: mainmaxsum.c
	$(CC) -c mainmaxsum.c $(CFLAGS)

maxsum.o: maxsum.c maxsum.h.gch
	$(CC) -c maxsum.c $(CFLAGS)

maxsum.h.gch: maxsum.h
	$(CC) -c maxsum.h

c_timer.o: c_timer.c
	$(CC) -c c_timer.c $(CFLAGS)

clean:
	rm -f *.out *.err *.o *.h.gch
