CXX=$(shell python3 -m korali.cxx --compiler)
KORALICXXFLAGS=$(shell python3 -m korali.cxx --cflags)
KORALILIBS=$(shell python3 -m korali.cxx --libs)

CXXFLAGS = $(KORALICXXFLAGS) -O3 -Wall -Wextra -Wno-unused-parameter

.PHONY: all clean

all: build/solver

build/solver: build/model.o build/korali_setup.o
	$(CXX) -o $@ $^ $(KORALILIBS)

build/%.o: src/%.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

clean:
	rm -f build/*.o build/solver
