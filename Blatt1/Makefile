all: A1 A2

%: %/%.cpp %/build
	g++ -O3 -ffast-math -march=native -I /usr/include/eigen3 $< -o %/$@
	cd %; ./$@
	python %.py

%/build:
	cd %
	mkdir -p build
clean:
	rm -f A1/A1
	rm -f A2/A2

.PHONY: clean
