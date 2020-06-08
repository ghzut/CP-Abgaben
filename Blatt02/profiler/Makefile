CXXFLAGS = -O3 -ffast-math -march=native

all: profilerExample.cpp Profiler.cpp Profiler.h
	g++ $(CXXFLAGS) profilerExample.cpp profiler.cpp -o profilerExample
	./profilerExample
	
# 'make clean' entfernt Zwischendateien
clean:
	rm -f ./profilerExample

# .PHONY gibt an, dass clean auch funktioniert, wenn eine Datei namens 'clean' existiert
.PHONY: clean