all: A1

A1: Aufgabe1_LU-Zerlegung/lu.cpp
	g++ -O3 -ffast-math -march=native -I /usr/include/eigen3 $< -o Aufgabe1_LU-Zerlegung/lu.out
	cd Aufgabe1_LU-Zerlegung; ./lu.out

clean:
	rm -f Aufgabe1_LU-Zerlegung/lu.out
	
.PHONY: clean