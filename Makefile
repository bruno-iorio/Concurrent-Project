CXX = g++                                # C++ compiler
CXXFLAGS = -std=c++20 -fopenmp -O3                  # Compiler flags

# Source and object files
SOURCES = main.cpp delta_step_dynamic.cpp delta_step_static.cpp dijkstra.cpp graph.cpp
OBJECTS = $(SOURCES:.cpp=.o)             # Creates main.o src/objects.o src/utils.o

EXEC = main                              # Output executable

# Default target
all: $(EXEC)

# Pattern rule to build object files
%.o: %.cpp
	$(CXX) -c $< -o $@

# Link the object files to create the executable
$(EXEC): $(OBJECTS)
	$(CXX)  $(OBJECTS) -o $(EXEC)

# Clean up build files
clean:
	rm -f $(OBJECTS) $(EXEC)

# Phony targets
