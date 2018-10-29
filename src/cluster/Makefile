TARGET ?= meshclust
VERSION ?= 1.1.0
CXX ?= g++
ifeq ($(debug),yes)
	CXXFLAGS += -ggdb  -fopenmp
else
	CXXFLAGS += -fopenmp -O3 -march=native -g
endif
CXXFLAGS += -std=c++11 -DVERSION=\"$(VERSION)\"
LDFLAGS += -lm

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
