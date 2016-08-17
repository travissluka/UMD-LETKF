.PHONY: build doc test clean all

all: build

clean:
	rm -rf build doc

build:
	mkdir -p build
	cd build; cmake ../ -DCMAKE_BUILD_TYPE=Debug; make

test: build
	cd build; ./letkf_test

doc:
	ford doc.md
