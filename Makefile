.PHONY: build

all: build

clean:
	rm -rf build

build:
	mkdir -p build
	cd build; cmake ../ -DCMAKE_BUILD_TYPE=Debug; make

test: build
	cd build; ./letkf_test
