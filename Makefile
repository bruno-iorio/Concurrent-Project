CXX = g++                             # C++ compiler
CXXFLAGS = -fopenmp -O3                  # Compiler flags
INCLUDES = -I.                           # Include current directory for .hpp files

# Source and object files
SOURCES = main.cpp  
OBJECTS = $(SOURCES:.cpp=.o)             # Creates main.o src/objects.o src/utils.o

EXEC = main                              # Output executable

# Default target
all: $(EXEC)

# Pattern rule to build object files
%.o: %.cpp
	$(CXX) $(INCLUDES) -c $< -o $@

# Link the object files to create the executable
$(EXEC): $(OBJECTS)
	$(CXX)  $(OBJECTS) -o $(EXEC)
# Clean up build files
clean:
	rm -f $(OBJECTS) $(EXEC)

# Phony targets
