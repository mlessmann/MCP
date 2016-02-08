CXXFLAGS := -std=c++11 -O2 -Wall -fopenmp
#CXXFLAGS += -pg        # Profiling von Funktionen: $ gprof benchmark
#CXXFLAGS += --coverage # Profiling von Codeblöcken: $ gcov benchmark.cpp

# Benchmark sources:
HEADERS := sequential.h parallel.h
SOURCES := benchmark.cpp

.PHONY: all
all: benchmark
	./$<

benchmark: $(HEADERS) $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: doc
doc: doc/Protokoll.tex
	$(MAKE) -C doc

.PHONY: clean
clean:
	$(MAKE) -C doc clean
	#rm -f benchmark
