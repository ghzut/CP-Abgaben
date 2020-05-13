A2: A2/A2.cpp A2/build 
	g++ -O3 -ffast-math -march=native -I /usr/include/eigen3 $< -o A2/$@
	cd A2; ./$@ > log.txt

A2/build:
	cd A2; mkdir -p build

clean:
	rm -f A2/A2
	rm -rf A2/build
	rm -f A2/log.txt