all: stable


stable: stable.cc
		g++ -O3 -ffast-math -march=native stable.cc -o stable
		./stable
		python stable.py

clean:
	  rm -f stable
