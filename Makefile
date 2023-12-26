CC = gcc
CFLAGS = -Wall -Wextra -g
DEPS = LLL_Reduction.h #Enumeration.h
SOURCE_FILES = LLL_Reduction.c main.c #Enumeration.c test.c
OBJECTS = LLL_Reduction.o main.o #Enumeration.o test.o

all: $(OBJECTS)
	$(CC) -o runme $(OBJECTS) -lm $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

test: $(OBJECTS)
	$(CC) -o testme $(OBJECTS) -lm $(CFLAGS)
	./testme

clean: 
	rm -rf runme testme *.o
