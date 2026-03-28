GTEST_PREFIX := $(shell brew --prefix googletest)

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Iinclude -Iexperiments
GTEST_INC = -I$(GTEST_PREFIX)/include
GTEST_LIB = -L$(GTEST_PREFIX)/lib -lgtest -lgtest_main -pthread

SRC = src/Matrix.cpp src/Gauss.cpp src/Generators.cpp experiments/Experiment.cpp
MAIN = src/main.cpp
TESTS = src/Tests.cpp

all: main

main:
	$(CXX) $(CXXFLAGS) $(MAIN) $(SRC) -o main

tests:
	$(CXX) $(CXXFLAGS) $(TESTS) $(SRC) $(GTEST_INC) $(GTEST_LIB) -o tests

clean:
	rm -f main tests