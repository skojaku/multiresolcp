# Makefile
.PHONY: all

#CC := gcc
CC := g++

#CFLAGS := -O3 -std=c++11 # use this option in case openmp does not work 
CFLAGS := -O3 -std=c++11 -fopenmp

all: km_ompnet

km_ompnet:src/* include/* km_ompnet/* 
	sudo rm -rf build km_multiresol.egg* && sudo python3 setup.py build install

.PHONY: clean
clean:
	$(RM) km_multiresol
