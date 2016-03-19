all: main

main: main.o
	g++ -std=c++11 -O3 main.o -o test -llapacke

main.o: main.cpp
	g++ -std=c++11 -O3 -c main.cpp -llapacke

clean:
	rm -rf *.o test *.gnu *.temp *.jpg load *~
