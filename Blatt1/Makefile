all: A1

A1: Aufgabe1_LU-Zerlegung/lu.cpp
	g++ -O3 -ffast-math -march=native -I /usr/include/eigen3 $< -o Aufgabe1_LU-Zerlegung/lu
	cd Aufgabe1_LU-Zerlegung; ./lu

clean:
	rm -f ./program
	
.PHONY: clean