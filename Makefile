CXX = g++
CXXFLAGS = -O3 -Wall -std=c++11
BIN = ABC++utils

all: $(BIN)

ABC++utils: ABC++utils.cpp ABC++utils.h
	$(CXX) $(CXXFLAGS) ABC++utils.cpp -o ABC++utils

clean:
	rm -f $(BIN) *.o *.d

.PHONY: clean all
