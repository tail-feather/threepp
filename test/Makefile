UNAME := $(shell uname)
UNAME_ARCH := $(shell uname -m)
ifeq ($(UNAME), Darwin)
	TARGET_OS := osx
else
	TARGET_OS := linux
endif
ifeq ($(TARGET_OS), osx)
	CXX := clang++
else
	CXX := g++-10
endif
BUILDLEVEL=-g
CXXFLAGS=-std=c++20 $(BUILDLEVEL) -Wall -fPIC -I./ -I../src/
SRCDIRS:=$(shell find * -type d)
SOURCES=$(foreach dir, $(SRCDIRS), $(wildcard $(dir)/*.cpp))
OBJECTS=$(addsuffix .o, $(basename $(SOURCES) main.cpp))

.PHONY: all

all: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o test $^
	#valgrind --leak-check=full ./test --log_level=test_suite
	./test --log_level=test_suite

clean:
	for file in $(OBJECTS); do if [ -f "$$file" ]; then rm "$$file"; fi; done
