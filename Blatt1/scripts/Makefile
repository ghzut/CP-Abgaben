all: two

two: two.cpp build
	g++ -I eigen3_3_7/ -O3 -ffast-math -march=native two.cpp -o two
	./two
	python two.py

build:
	mkdir build
