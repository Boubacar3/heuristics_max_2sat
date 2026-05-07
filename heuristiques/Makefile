# ==========================================
# Makefile for CPLEX C++ Project
# ==========================================

# 1. Compiler settings
CXX = g++

# 2. CPLEX Installation Directories
CPLEX_STUDIO_DIR = /opt/ibm/ILOG/CPLEX_Studio2212
CPLEXDIR         = $(CPLEX_STUDIO_DIR)/cplex
CONCERTDIR       = $(CPLEX_STUDIO_DIR)/concert

# 3. System and Architecture (Standard for 64-bit Linux)
SYSTEM    = x86-64_linux
LIBFORMAT = static_pic

# 4. Include directories (Headers)
INCLUDES = -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include

# 5. Compiler Flags
# -m64: Compile for 64-bit architecture
# -O3: Maximum optimization
# -fPIC: Position Independent Code (required by CPLEX)
# -fexceptions: Enable exception handling
# -DNDEBUG: Turn off debugging (for production speed)
# -DIL_STD: Tell CPLEX to use the standard C++ library
CXXFLAGS = -m64 -O3 -fPIC -fexceptions -DNDEBUG -DIL_STD $(INCLUDES)

# 6. Library directories
LIBDIRS = -L$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT) \
          -L$(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

# 7. Libraries to link against
LIBS = -lilocplex -lcplex -lconcert -lm -pthread -ldl

# 8. File gathering
# Find all .cpp files in the current folder
SRCS = $(wildcard *.cpp)
# Create a list of corresponding .o (object) files
OBJS = $(SRCS:.cpp=.o)

# 9. Executable name
TARGET = maxsat_solver

# ==========================================
# Build Rules
# ==========================================

# Default target: build the executable
all: $(TARGET)

# Rule to link the executable
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET)..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBDIRS) $(LIBS)
	@echo "Build successful! Run with ./$(TARGET)"

# Rule to compile .cpp files into .o object files
%.o: %.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to clean up build artifacts
clean:
	@echo "Cleaning up..."
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean