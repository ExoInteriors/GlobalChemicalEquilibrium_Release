# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++11

# Files
TARGET = solver
HEADERS = ../src/solver.h parameters.h
OBJECTS = ctot.o grad_ctot.o calculate_P.o optimizers.o solver.o main.o

all: solver

solver: $(OBJECTS) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o solver $(OBJECTS)

main.o: $(HEADERS) ../src/main.cpp
	$(CXX) $(CXXFLAGS) -c -o main.o ../src/main.cpp

solver.o: $(HEADERS) ../src/solver.cpp
	$(CXX) $(CXXFLAGS) -c -o solver.o ../src/solver.cpp

optimizers.o: $(HEADERS) ../src/optimizers.cpp
	$(CXX) $(CXXFLAGS) -c -o optimizers.o ../src/optimizers.cpp

ctot.o: $(HEADERS) ctot.cpp
	$(CXX) $(CXXFLAGS) -c -o ctot.o ctot.cpp

grad_ctot.o: $(HEADERS) grad_ctot.cpp
	$(CXX) $(CXXFLAGS) -c -o grad_ctot.o grad_ctot.cpp

calculate_P.o: $(HEADERS) calculate_P.cpp
	$(CXX) $(CXXFLAGS) -c -o calculate_P.o calculate_P.cpp


parameters.h ctot.cpp grad_ctot.cpp calculate_P.cpp : CTOT.py Equations.py
	python3 CTOT.py Equations.py

# Clean rule to remove generated files
clean:
	rm -f $(OBJECTS) $(TARGET)

# Convenience rule for rebuilding everything
rebuild: clean all

.PHONY: all clean rebuild
