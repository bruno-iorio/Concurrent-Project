CXX = g++                                # C++ compiler
CXXFLAGS = -std=c++20 -fopenmp -O3 -pthread  # Compiler flags

# Source and object files
SOURCES = main.cpp delta_step_dynamic.cpp delta_step_static.cpp dijkstra.cpp graph.cpp
OBJECTS = $(SOURCES:.cpp=.o)

EXEC = main                              # Output executable

# Default target
all: $(EXEC)

# Pattern rule to build object files with flags
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link the object files to create the executable with flags
$(EXEC): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXEC)

# Clean up build files
clean:
	rm -f $(OBJECTS) $(EXEC)

.PHONY: all clean
