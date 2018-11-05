CC=g++
CFLAGS=-c -g -O3 
SPLFLAG=-std=c++0x
all: execute

execute: main.o transpose.o parallel_multiply.o parallelGE.o parallelGS.o
	$(CC) -fopenmp  main.o transpose.o parallel_multiply.o parallelGE.o parallelGS.o -o execute

main.o: main.cpp
	$(CC) $(SPLFLAG) $(CFLAGS) main.cpp

transpose.o: transpose.cpp
	$(CC) $(CFLAGS) transpose.cpp

parallel_multiply.o: parallel_multiply.cpp
	$(CC) -fopenmp $(CFLAGS) parallel_multiply.cpp

gaussElimination.o: parallelGE.cpp
	$(CC) -fopenmp $(CFLAGS) parallelGE.cpp

gaussSiedel.o: parallelGS.cpp
	$(CC) -fopenmp $(CFLAGS) parallelGS.cpp

clean:
	rm *.o execute
