# Example commands:
#   make (build in release mode)
#   make debug (build in debug mode)
#   make clean (deletes *.o files, which aren't required to run the aligner)
#   make distclean (deletes *.o files and the binary)
#   make CXX=g++-5 (build with a particular compiler)
#   make CXXFLAGS="-Werror -g3" (build with particular compiler flags)


# CXX and CXXFLAGS can be overridden by the user.
CXX         ?= g++
CXXFLAGS    ?= -Wall -Wextra -pedantic -mtune=native

# These flags are required for the build to work.
LIB          = -lz
FLAGS        = -std=c++11

# Different debug/optimisation levels for debug/release builds.
DEBUGFLAGS   = -g
RELEASEFLAGS = -O3

TARGET       = bin/filtlong
SHELL        = /bin/sh
SOURCES      = $(shell find src -name "*.cpp")
HEADERS      = $(shell find src -name "*.h")
OBJECTS      = $(SOURCES:.cpp=.o)

.PHONY: release
release: FLAGS+=$(RELEASEFLAGS)
release: $(TARGET)

.PHONY: debug
debug: FLAGS+=$(DEBUGFLAGS)
debug: $(TARGET)

dir_guard=@mkdir -p $(@D)

$(TARGET): $(OBJECTS)
	$(dir_guard)
	$(CXX) $(FLAGS) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LIB)

clean:
	$(RM) $(OBJECTS)

distclean: clean
	$(RM) $(TARGET)

%.o: %.cpp $(HEADERS)
	$(CXX) $(FLAGS) $(CXXFLAGS) -c -o $@ $<
