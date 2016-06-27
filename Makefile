all: main

main: main.o
	g++ -g -std=c++11 -O3 main.o -o test -llapacke -lblas -fopenmp

main.o: main.cpp
	g++ -g -std=c++11 -O3 -c main.cpp -llapacke -lblas -fopenmp

clean:
	rm -rf *.o test *.gnu *.temp *.jpg load *.png *.dat *.vtk  *~
