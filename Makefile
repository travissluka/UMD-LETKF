.PHONY: build doc test clean all debug

all: build

clean:
	rm -rf build doc

build:
	mkdir -p build
	cd build; cmake ../src -DCMAKE_BUILD_TYPE=Release; make --no-print-directory

debug: 
	mkdir -p build
	cd build; cmake ../src -DCMAKE_BUILD_TYPE=Debug; make --no-print-directory

test: build
	cd build; make test

doc:
	ford doc.md
