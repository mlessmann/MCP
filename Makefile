CXXFLAGS := -std=c++11 -Wall -fopenmp

.PHONY: all
all: benchmark
	./$^

.PHONY: clean
clean:
	rm -f benchmark
