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
	$(CXX) $(CXXFLAGS) -c -o main.o ../src/main.cpp -I.

solver.o: $(HEADERS) ../src/solver.cpp
	$(CXX) $(CXXFLAGS) -c -o solver.o ../src/solver.cpp -I.

optimizers.o: $(HEADERS) ../src/optimizers.cpp
	$(CXX) $(CXXFLAGS) -c -o optimizers.o ../src/optimizers.cpp -I.

ctot.o: $(HEADERS) ctot.cpp
	$(CXX) $(CXXFLAGS) -c -o ctot.o ctot.cpp -I.

grad_ctot.o: $(HEADERS) grad_ctot.cpp
	$(CXX) $(CXXFLAGS) -c -o grad_ctot.o grad_ctot.cpp -I.

calculate_P.o: $(HEADERS) calculate_P.cpp
	$(CXX) $(CXXFLAGS) -c -o calculate_P.o calculate_P.cpp -I.


parameters.h ctot.cpp grad_ctot.cpp calculate_P.cpp : ../src/CTOT.py Equations.py ../Molecular_Weight.dat
	python3 ../src/CTOT.py Equations.py

# Clean rule to remove generated files
clean:
	rm -f $(OBJECTS) $(TARGET) parameters.h ctot.cpp grad_ctot.cpp calculate_P.cpp

# Convenience rule for rebuilding everything
rebuild: clean all

.PHONY: all clean rebuild
