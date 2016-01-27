CXXFLAGS := -std=c++11 -O2 -Wall -fopenmp
#CXXFLAGS += -pg        # Profiling von Funktionen: $ gprof benchmark
#CXXFLAGS += --coverage # Profiling von Codeblöcken: $ gcov benchmark.cpp
LDLIBS   := -lOpenCL

# Benchmark sources:
HEADERS := barrier.h sequential.h parallel.h
SOURCES := benchmark.cpp barrier.cpp

.PHONY: all
all: benchmark
	./$<

benchmark: $(HEADERS) $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: testcl
testcl: test-opencl
	./$<

.PHONY: clean
clean:
	rm -f benchmark
	rm -f test-opencl
