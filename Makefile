CC = gcc
CFLAGS = -Wall -Wextra -g
DEPS = LLL_Reduction.h Enumeration.h
SOURCE_FILES = LLL_Reduction.c Enumeration.c test.c main.c
OBJECTS = LLL_Reduction.o Enumeration.o test.o main.o

all: $(OBJECTS)
	$(CC) -o runme $(OBJECTS) -lm $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

test: $(OBJECTS)
	$(CC) -o testme $(OBJECTS) -lm $(CFLAGS)
	./testme

clean: 
	rm -rf runme testme *.o
