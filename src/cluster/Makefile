TARGET ?= meshclust
VERSION ?= 1.2.0

ifeq ($(shell uname), Darwin)
	CXX = xcrun -sdk macosx clang++
else
	CXX ?= g++
endif

ifeq ($(debug),yes)
	CXXFLAGS += -ggdb
else
	CXXFLAGS += -O3 -flto -march=native -DNDEBUG=1
endif
CXXFLAGS += -std=c++11 -pedantic -Wall -Wno-overloaded-virtual -DVERSION=\"$(VERSION)\"
LDFLAGS += -flto -lm

SOURCES := $(shell find ./src -name '*.cpp')
OBJECTS = $(SOURCES:%.cpp=bin/%.o)
BIN_OBJECTS := $(shell find  ../../bin/ -mindepth 2 -name '*.o')

all: clean $(TARGET)

$(TARGET): $(OBJECTS) $(BIN_OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

bin/%.o: %.cpp
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJECTS) $(TARGET)

install: $(TARGET)
	cp $(TARGET) ~/bin

.PHONY: all clean install
