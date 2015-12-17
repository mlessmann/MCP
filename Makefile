CXXFLAGS := -std=c++11 -Wall

.PHONY: all
all: sequential
	./$^

.PHONY: clean
clean:
	rm -f sequential
