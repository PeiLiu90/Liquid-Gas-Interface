DCF: main.o
	g++ -O3 -o DCF main.o
main.o: main.cpp FFT.h FMT2.h Gauss.h GMRES.h LinearEquation.h NewtonGmres.h parameter.h RHNC.h vector2d.h vector3d.h input.dat
	g++ -O3 -c main.cpp -o main.o

clean:
	rm DCF main.o
