all: A2 #A3 # A1


A1: A1/A1.cpp A1/build A1/A1.py
	g++ -O3 -ffast-math -march=native -I /usr/include/eigen3 $< -o A1/$@
	cd A1; ./$@ > log.txt
	python A1/A1.py


A2: A2/A2.cpp A2/build A2/A2.py
	g++ -O3 -ffast-math -march=native -I /usr/include/eigen3 -Wall $< -o A2/$@
	cd A2; ./$@ > log.txt
	python A2/A2.py

A1/build:
	cd A1; mkdir -p build


A2/build:
	cd A2; mkdir -p build

clean:
	rm -f A1/A1
	rm -f A2/A2
	rm -rf A1/build
	rm -rf A2/build
	rm -f A1/log.txt
	rm -f A2/log.txt
.PHONY: clean
