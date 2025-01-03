CC = gcc
CFLAGS = -Wall -Wextra -g
DEPS = LLL_Reduction.h Enumeration.h Threads.h
SOURCE_FILES = LLL_Reduction.c Enumeration.c GeneralFunctions.c main.c test.c
OBJECTS = LLL_Reduction.o Enumeration.o GeneralFunctions.o

all: $(OBJECTS) main.o
	$(CC) -o runme $(OBJECTS) main.o -lm $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

test: test.o $(OBJECTS)
	$(CC) -o testme $(OBJECTS) test.o -lm $(CFLAGS)
	./testme

clean: 
	rm -rf runme testme *.o result.txt
