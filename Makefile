CXXFLAGS := -std=c++11 -Wall \
            -fopenmp
#CXXFLAGS += -pg        # Profiling von Funktionen: $ gprof benchmark
#CXXFLAGS += --coverage # Profiling von Codeblöcken: $ gcov benchmark.cpp
LDLIBS   := -lOpenCL

.PHONY: all
all: benchmark
	./$<

.PHONY: testcl
testcl: test-opencl
	./$<

.PHONY: clean
clean:
	rm -f benchmark
	rm -f test-opencl
