NAME ?= meanshift
CXX ?= g++
ifeq ($(debug),yes)
	CXXFLAGS += -ggdb
else
	CXXFLAGS += -O3 -march=native
endif
CXXFLAGS += -fopenmp -std=c++11
LDFLAGS += -lm

SOURCES := $(shell find ./src -name '*.cpp')
OBJECTS = $(SOURCES:%.cpp=bin/%.o)
BIN_OBJECTS := $(shell find  ../../bin/ -mindepth 2 -name '*.o')

all: clean $(NAME)

$(NAME): $(OBJECTS) $(BIN_OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

bin/%.o: %.cpp
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJECTS) $(NAME)

install: $(NAME)
	cp $(NAME) ~/bin

.PHONY: all clean install
