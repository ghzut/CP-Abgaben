all: 0_1a.pdf 0_1b.pdf 0_1c.pdf

%.pdf: stable.cc
		g++ -O3 -ffast-math -march=native stable.cc -o stable
		./stable
		python stable.py
		
clear:
		rm -rf stable
