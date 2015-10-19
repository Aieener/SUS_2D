# declare the variable
# CC = g++ -std=gnu++0x
CC = g++

CFLAGS = -c -Wall -pedantic -Ofast -march=native -std=c++11 

PROG1 = main

all: runit

runit: $(PROG1).o cells.o square.o hardrods.o histogram.o MC.o 
		$(CC) $(PROG1).o cells.o square.o hardrods.o histogram.o MC.o -o runit

MC.o: MC.cpp
		$(CC) $(CFLAGS) MC.cpp

$(PROG1).o: $(PROG1).cpp
		$(CC) $(CFLAGS) $(PROG1).cpp

cells.o: cells.cpp
		$(CC) $(CFLAGS) cells.cpp

square.o: square.cpp
		$(CC) $(CFLAGS) square.cpp

hardrods.o: hardrods.cpp
		$(CC) $(CFLAGS) hardrods.cpp

histogram.o: histogram.cpp
		$(CC) $(CFLAGS) histogram.cpp

clean:
		rm -rf *o run

