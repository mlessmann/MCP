CXXFLAGS := -std=c++11 -Wall \
            -fopenmp
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
