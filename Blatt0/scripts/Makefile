all: stable

stable: stable.cc build
		g++ -O3 -ffast-math -march=native stable.cc -o build/stable
		./build/stable
		python stable.py

build:
	mkdir -p build

clean:
	rm -rf build
