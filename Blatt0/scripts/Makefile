ifeq (,$(shell sh -c 'cygpath --version 2> /dev/null'))
  # Unix
  pwd := $$(pwd)
  translate = $1
else
  # Windows mit MSys2/Cygwin
  pwd := $$(cygpath -m "$$(pwd)")
  translate = $(shell echo '$1' | sed 's/:/;/g')
endif

all: stable

stable: stable.cc build
		g++ -O3 -ffast-math -march=native stable.cc -o build/stable
		./build/stable
		python stable.py

build:
	mkdir -p build

clean:
	rm -rf build

