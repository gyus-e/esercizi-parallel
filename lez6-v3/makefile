a.out: main.o maxsum.o
	gcc -o a.out main.o maxsum.o -fopenmp -lm

main.o: main.c maxsum.h
	gcc -c main.c -fopenmp

maxsum.o: maxsum.c maxsum.h
	gcc -c maxsum.c maxsum.h -fopenmp

clean:
	rm a.out *.o
