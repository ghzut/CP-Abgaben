all: A1 A2 A3


A1: A1/A1.cpp A1/build A1/A1.py
	g++ -O3 -ffast-math -march=native -I /usr/include/eigen3 $< -o A1/$@
	cd A1; ./$@ > log.txt
	python A1/A1.py


A2: A2/A2.cpp A2/build
	g++ -O3 -ffast-math -march=native -I /usr/include/eigen3 -fopenmp $< -o A2/$@
	cd A2; ./$@ > log.txt
	python A2/A2.py

A3: A3/A3.cpp A3/build A3/A3.py
	g++ -O3 -ffast-math -march=native -I /usr/include/eigen3 -fopenmp $< -o A3/$@
	cd A3; ./$@ > log.txt
	python A3/A3.py


A1/build:
	cd A1; mkdir -p build


A2/build:
	cd A2; mkdir -p build

A3/build:
	cd A3; mkdir -p build

clean:
	rm -f A1/A1
	rm -f A2/A2
	rm -f A3/A3
	rm -rf A1/build
	rm -rf A2/build
	rm -rf A3/build
	rm -f A1/log.txt
	rm -f A2/log.txt
	rm -f A3/log.txt
.PHONY: clean
