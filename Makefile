CXXFLAGS := -std=c++11 -Wall

.PHONY: all
all: benchmark
	./$^

.PHONY: clean
clean:
	rm -f benchmark
