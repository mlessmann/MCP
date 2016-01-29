CXXFLAGS := -std=c++11 -O2 -Wall -fopenmp
#CXXFLAGS += -pg        # Profiling von Funktionen: $ gprof benchmark
#CXXFLAGS += --coverage # Profiling von Codeblöcken: $ gcov benchmark.cpp
OCLLIBS  := -lOpenCL

# Benchmark sources:
HEADERS := sequential.h parallel.h
SOURCES := benchmark.cpp

.PHONY: all
all: benchmark
	./$<

benchmark: $(HEADERS) $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: test-opencl
test-opencl: opencl
	./$<

opencl: opencl.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(OCLLIBS)

.PHONY: clean
clean:
	rm -f benchmark
	rm -f opencl
