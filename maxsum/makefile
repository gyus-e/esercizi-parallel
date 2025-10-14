CC = gcc
CFLAGS = -O3 -fopenmp -Wall -Wextra -std=gnu90
LDFLAGS = -O3 -fopenmp
LDLIBS = -lm

mainmaxsum_v2.out: mainmaxsum.o maxsum_v2.o c_timer.o mainmaxsum.out
	$(CC) -o mainmaxsum_v2.out mainmaxsum.o maxsum_v2.o c_timer.o $(LDFLAGS) $(LDLIBS)

mainmaxsum.out: mainmaxsum.o maxsum.o c_timer.o
	$(CC) -o mainmaxsum.out mainmaxsum.o maxsum.o c_timer.o $(LDFLAGS) $(LDLIBS)

mainmaxsum.o: mainmaxsum.c
	$(CC) -c mainmaxsum.c $(CFLAGS)

maxsum_v2.o: maxsum_v2.c
	$(CC) -c maxsum_v2.c $(CFLAGS)

maxsum.o: maxsum.c
	$(CC) -c maxsum.c $(CFLAGS)

c_timer.o: c_timer.c
	$(CC) -c c_timer.c $(CFLAGS)

clean:
	rm -f *.out *.err *.o
