CXXFLAGS := -std=c++11 -Wall -fopenmp -O2

.PHONY: all
all: benchmark
	./$^

.PHONY: clean
clean:
	rm -f benchmark
