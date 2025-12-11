CC = gcc
OPTIONS = -O3 -Wall -Wextra -std=gnu90
CFLAGS = -fopenmp 
LDFLAGS = -fopenmp
LDLIBS = -lm

mainmaxsum_v2.out: mainmaxsum.o maxsum_v2.o c_timer.o mainmaxsum.out
	$(CC) $(OPTIONS) -o mainmaxsum_v2.out mainmaxsum.o maxsum_v2.o c_timer.o $(LDFLAGS) $(LDLIBS)

mainmaxsum.out: mainmaxsum.o maxsum.o c_timer.o
	$(CC) $(OPTIONS) -o mainmaxsum.out mainmaxsum.o maxsum.o c_timer.o $(LDFLAGS) $(LDLIBS)

mainmaxsum.o: mainmaxsum.c
	$(CC) $(OPTIONS) -c mainmaxsum.c $(CFLAGS)

maxsum_v2.o: maxsum_v2.c
	$(CC) $(OPTIONS) -c maxsum_v2.c $(CFLAGS)

maxsum.o: maxsum.c
	$(CC) $(OPTIONS) -c maxsum.c $(CFLAGS)

c_timer.o: c_timer.c
	$(CC) $(OPTIONS) -c c_timer.c $(CFLAGS)

clean:
	rm -f *.out *.err *.o
