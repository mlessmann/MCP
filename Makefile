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
	# Mehrgitter Plot:
	cp Mehrgitter.plot.template Mehrgitter.plot
	for file in *.matrix; do \
		echo "set output '$$file.png'" >> Mehrgitter.plot; \
	    echo "splot '$$file' matrix using 1:2:3 with lines" >> Mehrgitter.plot; \
	done
	gnuplot Mehrgitter.plot

benchmark: $(HEADERS) $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: test-opencl
test-opencl: opencl
	./$<

opencl: opencl.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(OCLLIBS)

.PHONY: doc
doc: doc/Protokoll.tex
	$(MAKE) -C doc

.PHONY: clean
clean:
	$(MAKE) -C doc clean
	#rm -f benchmark
	#rm -f opencl
	rm -f *.matrix
	rm -f Mehrgitter.plot
